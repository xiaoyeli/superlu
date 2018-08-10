find_package(Matlab)
# This module is untested because I don't have Matlab
# but I want to provide something for completeness and
#hopefully, someone else will work on this..
if(Matlab_FOUND)
    matlab_add_mex(mexlusolve SRC mexlusolve.c  LINK_TO superlu )
    matlab_add_mex(mexlusolve SRC mexsuperlu.c  LINK_TO superlu )
else()
    message(STATUS WARNING "This requires MATLAB from Mathworks.")
endif(Matlab_FOUND)
