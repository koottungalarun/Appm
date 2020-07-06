# Run the APPM application
#
# If the PowerShell script cannot be run: check your execution policy 
# using the cmdlet
#
# Get-ExecutionPolicy
# 
# See also: https://docs.microsoft.com/de-de/powershell/module/microsoft.powershell.core/about/about_execution_policies?view=powershell-7
# 
# You might have to edit the execution policy:
#
# Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
#
#
# The executable requires that Intel MKL is available on the 
# environment path. If an error is shown saying that 'mkl_sequential.dll' 
# and/or 'mkl_core.dll' is not found, add Intel MKL to your environment path.
#
# Example: 
# $env:Path = $env:Path + ';C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\intel64\mkl'
#
#
# 


####################################################
# Define application
$app = "..\x64\Release\appm.exe"

####################################################
# Run application (&-command), 
# and redirect output in two directions
#  (i) to console, and 
# (ii) to log file 
#

& $app | Tee-Object -FilePath "out.log" 

#
#
#
######################
# TODO
######################
# Define a database folder 
# copy executable to that folder
# create subfolders for testcases
# copy required input files to each testcases
# run testcases
# evaluate results with Paraview and save it with state file
# when loading the statefile into paraview, select correctly how to read
# adapt settings and save data
