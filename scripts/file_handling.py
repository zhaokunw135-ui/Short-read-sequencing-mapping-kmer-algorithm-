class FileHandling:
    def __init__(self, file_name):
        """
        Arg:
            file_name: transferred to create attribute for self (instance)
        """
        self.file_name = file_name

        self.check_file_exist()
        self.check_file_format()
        self.check_file_sequence()

    def check_file_exist(self):
        """
        check whether the file exists
        """
        try:
            with open(self.file_name, 'r'):
                print(f"File {self.file_name} found.")
        except FileNotFoundError:
            raise FileNotFoundError(f"Error! Can NOT find file'{self.file_name}'.")

    def check_file_format(self):
        """
        check whether the file is valid (first line starts with ">")
        """
        with open(self.file_name, 'r') as file:
            first_line = file.readline().strip()
            if not first_line.startswith(">"):
                raise ValueError(f"Error! '{self.file_name}' is invalid in format")

    def check_file_sequence(self):
        """
        check whether the sequences are valid (contain invalid characters)
        """
        with open(self.file_name, 'r') as file:
            valid_chars = set("AaTtGgCcNn")
            for line in file:
                line = line.strip()
                if not line.startswith(">"):
                    if not set(line).issubset(valid_chars):
                        raise ValueError(f"Error! In '{self.file_name}' line  {line} contains invalid character")
                    continue

    def load_fasta(self):
        """
            Reads a FASTA file and returns a dictionary of sequences.
            Return:
                load: A dictionary in which keys are sequence headers (remove '>')
            and values are the corresponding sequences.
        """
        load = {}
        with open(self.file_name, 'r') as file:
            header = None
            for line in file:
                line = line.strip()
                if line.startswith(">"):  # recognized the header
                    if header:
                        load[header] = "".join(sequence).upper()  # store the set before
                    header = line[1:]  # initialize a new set
                    sequence = []
                else:  # recognized the sequence
                    sequence.append(line)
            if header:
                load[header] = "".join(sequence).upper()  # store the last set
        print(f"File {self.file_name} is loaded.")
        return load
