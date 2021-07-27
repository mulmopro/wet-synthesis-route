'''----------------------------------------------------------------------------
Copyright © 2021 Politecnico di Torino

This file is part of WetSynthRoute.

WetSynthRoute is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WetSynthRoute is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WetSynthRoute.  If not, see <https://www.gnu.org/licenses/>.
----------------------------------------------------------------------------'''


import os
import sys
# import numpy as np


# Class derived from Exception for general invalid input in "caseSetup" file
class InputDictErr(Exception):

    def __init__(self, dictName, errorType, lineNum=None, line=None,
                 keyName=None):
        if errorType == 'dictNotFound':
            self.message = 'Dictionary "' + dictName + \
                            '" is not defined in the "caseSetup" file.\n' + \
                            'Please specify required data for the "' + \
                            dictName + '".\n'

        elif errorType == 'inputNotFound':
            self.message = '"' + dictName + \
                            '" is not defined in the "caseSetup" file.\n' + \
                            'Please specify required data for the "' + \
                            dictName + '".\n'

        elif errorType == 'keyError':
            self.message = '"' + keyName + '" is not defined in the "' + \
                            dictName + '" dictionary.\n' \
                            'Please check the input file.\n'

        elif errorType == 'invalidData':
            self.message = 'Data specified for "' + keyName + '" in "' + \
                            dictName + '" dictionary is not valid:\n' + \
                            line + '\nPlease check the input file.\n'

        elif errorType == 'invalidModel':
            self.message = 'Invalid "' + keyName + \
                            '" model is specified in "' + dictName + \
                            '" dictionary:\n' + line + \
                            '\nPlease check the input file.\n'

        elif errorType == 'invalidType':
            self.message = 'Invalid "' + keyName + \
                            '" is specified in "' + dictName + \
                            '" dictionary:\n' + line + \
                            '\nPlease check the input file.\n'

    def __str__(self):
        return '\n\nError while reading input files:\n' + self.message


# Class derived from Exception for realizability issue
class RealizabilityErr(Exception):

    def __init__(self, moments, negValue=None, negNode=False):
        if negNode:
            self.message = '\nError:\nNegative node "' + str(negValue) + \
                '" is calculated for the following set of moments:\n' + \
                str(moments)

        else:
            self.message = '\nError:\nRealizability issue with the ' + \
                'following set of moments:\n' + str(moments)

    def __str__(self):
        return self.message


# Class derived from FileNotFoundError to be raised
# in the case of not finding the required file
class FileNotFoundErr(FileNotFoundError):

    def __init__(self, filePath):
        dirName, fileName = os.path.split(filePath)
        self.message = '\nThe file "' + fileName + \
            '" does not exist in location:\n' + \
            '"' + dirName + '". Please check.\n'

    def __str__(self):
        return self.message


# Class derived from Exception to stop the ode solver
# if the ode function is called more than maxFuncEval
class MaxFuncEvalErr(Exception):

    def __init__(self, maxFuncEval):
        self.message = '\nError:\nThe ode solver is stopped after ' + \
            str(maxFuncEval) + ' number of function evaluations.\n'

    def __str__(self):
        return self.message


# Class derived from Exception to be raised when the integration fails
class IntegrationFailedErr(Exception):

    def __init__(self):
        self.message = '\nError:\nIntegration step failed.\n'

    def __str__(self):
        return self.message


# Function to modify string property of unforeseen exceptions raised in runtime
def create_unhandled_exception(raisedException):

    return type(raisedException)(
        '\nUnhandled exception happened: \n' + str(raisedException) + '\n'
        ).with_traceback(sys.exc_info()[2])
