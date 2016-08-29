import subprocess
import os


class SimpleExecutable(object):
    """Wrapper for external executable.

    .. py::attribute:: name
        (string) The name of the program (human readable).

    .. py::attribute:: path
        (string) The path pointing to the executable.

    .. py::attribute:: description
        (None or string) Describing the function of the executable,
        possibly with input and output specified.

    .. py::attribute:: parameters
        (None or function) Should be a function taking one argument,
        returning a list of strings. This list is then passed as
        command-line arguments to the executable. This function should
        be able to handle `None` as an argument.
    """
    def __init__(self, name, description, path,
                 working_dir=None, parameters=None):
        self.name = name
        self.path = path
        self.description = description
        self.working_dir = working_dir
        self.parameters = parameters

    def run(self, args_obj=None, **kwargs):
        """Call `subprocess.run`.

        :param args_obj:
            Object containing information for arguments. This is passed
            through the :py:attribute:`parameters` function attribute to
            generate the list of command-line arguments.
        :type args_obj: Any

        :param **kwargs:
            Keyword arguments are passed to `subprocess.run`.

        :return:
            CompletedProcess object.
        """
        if self.working_dir:
            orig_wd = os.getcwd()
            os.chdir(self.working_dir)

        args = [self.path]
        if self.parameters:
            args.extend(self.parameters(args_obj))

        result = subprocess.run(args, **kwargs)

        if self.working_dir:
            os.chdir(orig_wd)

        return result

