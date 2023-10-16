import os

class OutputFile:
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.file_name = os.path.splitext(os.path.basename(file_path))[0]

    def write_head(self, column_list: list):
        head_str = ''
        for column in column_list:
            head_str += column + '\t'
        with open(self.file_path, 'w') as file:
            file.write(head_str.rstrip() + '\n')

    def write_line(self, line: str):
        with open(self.file_path, 'a') as file:
            file.write(line)