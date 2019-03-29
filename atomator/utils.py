import os
import os.path
import re

def get_file_extension(filename):
    """Return the extension of a file.

        Arguments:
         - filename - string, the file to extract the extension
    
    """
    # If we want to get compound extensions, we just need to add it in the reg exp
    return re.compile(r'^.*?[.](?P<ext>fa\.gz|fasta\.gz|vcf\.gz|\w+)$').match(filename).group('ext')

def get_files(input_path, admited_formats):
    """Return a list of path files.

        Arguments:
         - input_path - string, the path where will look for the files
         - admited_formats - set, the admited formats for the files

    """
    result = []
    folder = "."
    if os.path.isdir(input_path):
        folder = input_path
        result = [os.path.abspath(os.path.join(folder, f)) for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f)) and get_file_extension(os.path.join(folder,f)) in admited_formats ]
        result.sort()
    else:
        filename = os.path.join(folder, input_path)
        if os.path.isfile(filename) and get_file_extension(filename) in admited_formats:
            result = [os.path.abspath(filename)]
    return result