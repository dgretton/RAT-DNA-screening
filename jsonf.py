import json
import os

def jsonf(fname):
    with open(fname) as f:
        return json.loads(f.read())

def jsonw(struct, fname, indent=2):
    with open(fname, 'w+') as f:
        f.write(json.dumps(struct, indent=indent))

def json_is_equivalent(json1, json2):
    return json.dumps(json1, sort_keys=True) == json.dumps(json2, sort_keys=True)

def json_deep_copy(json_in):
    return json.loads(json.dumps(json_in))

def json_in(json_sub, json_super):
    if json_is_equivalent(json_super, json_sub):
        return True
    if type(json_sub) is list:
        if type(json_super) is not list:
            return False
        if len(json_sub) > len(json_super):
            return False
        for super_i, sub_i in zip(json_super, json_sub):
            if super_i == sub_i:
                continue
            if json_in(super_i, sub_i):
                continue
            return False
        return True
    elif type(json_sub) is dict:
        if type(json_super) is not dict:
            return False
        for k, v in json_sub.items():
            if k not in json_super:
                return False
            if json_super[k] == v:
                continue
            if json_in(json_super[k], json_sub[k]):
                continue
            return False
        return True
    else:
        # it was a value, but not json-equivalent
        return False

# the meta data dictionary layers:
# level one, keys: sequences, sequence_meta, meta.
# sequences: value = a list of all the sequences.
# sequences_meta: keys: each sequence, value: a list, which contains a dict:
#### the keys : hazard, index, reverse complement, and so on.
# within meta, the keys are: creation time and so on.

# step 1: initialize a dict with the above characteristics
# step 2: add a sequence: this would have to add to both the sequences list and the sequence_meta, with the index it's added at in the sequence keys
# step 3: add meta data for a specified sequence, given the sequence and the meta data descriptor
# step 4: add general meta data, included under the meta key in the dict. check with dana about this.

def json_initialize_meta_dict():
    # initializes backbone structure, ie, the meta_dict with the three main keys.
    meta_dict= {}
    meta_dict['sequences']=[]
    meta_dict['sequence_meta']={}
    meta_dict['meta']={}
    return meta_dict

def write_empty_fragset(fname, overwrite=False):
    if not overwrite and os.path.isfile(fname):
        raise FileExistsError(f'{fname} already exists. pass overwrite=True to force.')
    jsonw(json_initialize_meta_dict(), fname)

def json_add_sequence(meta_dict, sequence):
    # adds a given sequence to the meta_dict. adds the sequence to the sequences value list, sequence as a key t
    # to the sequence_meta key value dict as a key, with the value that is a dict with the index as one of the keys.

    index = len(meta_dict['sequences']) + 1
    meta_dict['sequences'].append(sequence)

    meta_dict['sequence_meta'][sequence]= []
    sequence_meta_dict={}
    sequence_meta_dict["index"]=index

    meta_dict['sequence_meta'][sequence].append(sequence_meta_dict)

def json_add_meta_sequence(meta_dict, sequence, seq_meta):
    #adds meta data given the meta_type, the value of the meta data, and the specified sequence.

    if sequence in meta_dict['sequence_meta']:
        existing_metas = meta_dict['sequence_meta'][sequence]
        for existing_meta in existing_metas:
            if json_in(seq_meta, existing_meta):
                return # this metadata is fully & recursively already included
        existing_metas.append(json_deep_copy(seq_meta))
    else:
        json_add_sequence(meta_dict, sequence)
        meta_dict['sequence_meta'][sequence][0].update(json_deep_copy(seq_meta))

def json_add_meta_general(meta_dict, meta_type, meta_value):
    meta_dict['meta'][meta_type] = meta_value

def json_merge_into(meta_dict, sequences=None, sequence_meta=None, meta=None):
    if sequences:
        for var_seq in sequences:
            json_add_sequence(meta_dict, var_seq)
    if sequence_meta:
        for var_seq, seq_meta_list in sequence_meta.items():
            for seq_meta in seq_meta_list:
                json_add_meta_sequence(meta_dict, var_seq, seq_meta)
    if meta:
        for key, val in meta.items():
            json_add_meta_general(meta_dict, key, val)

def json_merge_fragset_into(new_fragset, fragset):
        jsonf.json_merge_into(fragset, sequences=new_fragset['sequences'],
                sequence_meta=new_fragset['sequence_meta'], meta=new_fragset['meta'])

def json_get_all_meta(meta_dict, sequence):
    return meta_dict['sequence_meta'][sequence]

def json_get_specific_meta(meta_dict, sequence, meta_type):
    individual_sequence_meta_dict = ['sequence_meta'][sequence][0]
    return individual_sequence_meta_dict[meta_type]

def json_copy_dict(fname_copy_from, fname_copy_to):
    # copies the meta_dict in one file to the meta_dict in a new file so no original information is lost. 
    with open(fname_copy_from) as f, open(fname_copy_to) as to:
        to.write(f.read())


if __name__ == '__main__':
    import sys
    if '--fix' in sys.argv:
        fname = sys.argv[sys.argv.index('--fix') + 1]
        if fname.split('.')[-1] != 'json' and '-f' not in sys.argv:
            print('use -f to fix file without json extension', fname)
            exit()
        jsonw(jsonf(fname), fname)
