import os, sys, shutil

class BaseClass():

    def _error(self, message):
        """
            Method for printing errors in nice manner.
        """

        msg = "\n+" + "-" * 78 + "+" + "\n" + "| Datgen Error: "
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (82 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
 
        print(msg, flush=True)

        exit()

    def _initialization(self, type, options):
        """
            Standard method for performing initialization.
        """
        
        if options is not None:
            if not isinstance(options, dict):
                self._error("The 'options' argument provided is not a dictionary.")
            elif options == {}:
                self._error("The 'options' argument provided is an empty dictionary.")
        else:
            self._error("Options argument not provided.")

        if type == "multi":
            self._setupMultiAnalysis(options)
        elif type == "single":
            self._setupSingleAnalysis(options)

        else:
            self._error("Value of type argument not recognized.")

    def _getDefaultOptions(self, defaultOptions):
        """
            Setting up the initial values of options.
        """

        for key in vars(defaultOptions):
            value = getattr(defaultOptions, key)
            self.options[key] = value

    def _setOptions(self, options):
        """
            Method for assigning user provided options.
        """

        for key in options.keys():
            if isinstance(options[key], dict):
                if key in self.options.keys():
                    # If the value is dictionary, update the default dictionary.
                    # Otherwise, assign values.
                    self.options[key].update(options[key]) 
                else:
                    self.options[key] = options[key]
            else:
                self.options[key] = options[key]

    def _setDirectory(self):
        """
            Method for setting up directory.
        """

        directory = self.options["directory"]

        if not os.path.isdir(directory):
            os.system("mkdir {}".format(directory))
        else:
            os.system("rm -r {}".format(directory))
            os.system("mkdir {}".format(directory))

        if self.options["type"] == "multi":
            FFD = False
            for key in self.options["varyingParameters"].keys():
                if key == "twist" or key == "shape":
                    FFD = True
                    break

            for sampleNo in range(self.options["numberOfSamples"]):
                os.system("mkdir {}/{}".format(directory,sampleNo))
                pkgdir = sys.modules["datgen"].__path__[0]
                filepath = os.path.join(pkgdir, "runscripts/runscript_aerodynamics.py")
                shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
                os.system("cp -r {} {}/{}/grid.cgns".format(self.options["aeroSolverOptions"]["gridFile"],directory,sampleNo))

                # copying ffd file
                if FFD:
                    pkgdir = sys.modules["datgen"].__path__[0]
                    filepath = os.path.join(pkgdir, "runscripts/deform_mesh.py")
                    shutil.copy(filepath, "{}/{}".format(directory,sampleNo))
                    os.system("cp -r {} {}/{}/ffd.xyz".format(self.options["ffdFile"],directory,sampleNo))

