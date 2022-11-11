import os

def mkdir(dir_):
    os.makedirs(dir_, exist_ok=True)

def list_files(dir_, pattern = None, full_name = True):
    ret = []
    if pattern != None: 
        ret = [f for f in os.listdir(dir_) if f.find(pattern) != -1]
    else:
        ret = [f for f in os.listdir(dir_)]

    if full_name:
        ret = [dir_ + "/" + f for f in ret]

    return ret
