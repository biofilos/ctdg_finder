import os

def create_folders(folder, folders):
    """
    Create the directories where all the
    results will be saved
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
        for sub_dir in folders:
            target = folder + '/' + sub_dir
            if not os.path.exists(target):
                os.makedirs(target)
    print("Directory tree for {} was created.".format(folder))
