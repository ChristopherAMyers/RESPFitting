from distutils.util import strtobool

class _KeywordOption():
    def __init__(self, key, var, type=str) -> None:
        pass

class Options():
    def __init__(self, input_file=None) -> None:
        """ Read input file options

            Parameters
            ----------
            input_file: str
                location of the input file
        """
        self.name = ""
        self.basis = "6-31G*"
        self.method = "HF"
        self.optimize_first = False
        self.charge = 0
        self.amber_fitting = True
        self.lone_pairs = False

        self.density_fitting = True
        self.lone_pairs = False
        self.n_dens = 1
        self.nh_dens = 1

        #   keywords that need to be converted to types other than strings
        self._keywords = {}
        self._keywords['charge'         ] = int
        self._keywords['optimize'       ] = strtobool
        self._keywords['density_fitting'] = strtobool
        self._keywords['lone_pairs'     ] = strtobool

        self.input_sections = {}
        if input_file is not None:
            self.input_sections = self._read(input_file)
        else:
            print("NONE")
            #   no file provided means use the default options
            pass
    
    def _read(self, input_file):
        input_sections = {}
        reading_sec = None
        for line in open(input_file, 'r'):     
            line = line.strip()
            if "$" in line:
                if "$end" in line:
                    reading_sec = None
                else:
                    reading_sec = str.lower(line[1:])
                    input_sections[reading_sec] = []
            
            elif reading_sec is not None:
                sp_comment = line.split('!')
                if sp_comment[0] != '!':
                    sp = line.replace('=', '').split()
                    if len(sp) < 2:
                        raise ValueError("Invalid line in input file\n %s" % line)
                    input_sections[reading_sec].append(sp)

        #   dynamically assign and convert options if needed
        for line in input_sections['resp']:
            option = line[0].lower()
            value = line[1].lower()
            converter = self._keywords.get(option, str)
            setattr(self, option, converter(value))

        return input_sections



