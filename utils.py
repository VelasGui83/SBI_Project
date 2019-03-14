import os
import os.path
import re

def get_file_extension(filename):
    """
    
    """
    return re.compile(r'^.*?[.](?P<ext>fa\.gz|fasta\.gz|vcf\.gz|\w+)$').match(filename).group('ext')

def get_files(input_path, admited_formats):
    """

    """
    result = []
    folder = "."
    if os.path.isdir(input_path):
        folder = input_path
        result = [os.path.join(folder, f) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and get_file_extension(os.path.join(folder,f)) in admited_formats ]
        result.sort()
    else:
        filename = os.path.join(folder, input_path)
        if os.path.isfile(filename) and get_file_extension(filename) in admited_formats:
            result = [filename]
    return result

def remove_items_containing(string_list, keyword):
    """

    """
    string_list_copy = string_list.copy()
    for item in string_list_copy:
        if keyword in item:
            string_list.remove(item)