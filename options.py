from distutils.util import strtobool
from tokenize import String
import sys

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
        self.optimize_first = False
        self.charge = 0
        self.amber_fitting = True
        self.lone_pairs = False
        self.vdw_ratios = [1.4, 1.6, 1.8, 2.0]
        self.mk_density = 10

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
        self._keywords['n_dens'         ] = int
        self._keywords['nh_dens'        ] = int
        self._keywords['vdw_ratios'     ] = self.strToFloatList
        self._keywords['mk_density'     ] = int

        self.input_sections = {}
        if input_file is not None:
            self.input_sections = self._read(input_file)
        else:
            print("NONE")
            #   no file provided means use the default options
            pass


    def strToFloatList(self, string):
        for char in ['[', ']', '(', ')', '{', '}']:
            string = str.replace(string, char, '')
        string = str.replace(string, ',', ' ')
        return [float(x) for x in string.split()]

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
                    input_sections[reading_sec] = {}
            
            elif reading_sec is not None:
                #   skip over comment lines
                sp_comment = line.split('!')
                if sp_comment[0] == '!': continue
                #   comments can also be placed at the end of lines
                if len(sp_comment) == 2:
                    line = sp_comment[0].rstrip()

                #   split only to check the number of entries in the line
                sp = line.replace('=', '').split()
                if len(sp) < 2:
                    raise ValueError("Invalid line in input file\n %s" % line)
                option = sp[0]
                value = line.replace(sp[0], '').lstrip()
                input_sections[reading_sec][option] =  value

        #   dynamically assign options and convert types if needed
        for option, value in input_sections['resp'].items():
            option = option.lower()
            value = value.lower()
            converter = self._keywords.get(option, str)
            setattr(self, option, converter(value))

        return input_sections



