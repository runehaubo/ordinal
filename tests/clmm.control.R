library(ordinal)
data(wine)


### 3 options for specifying control arguments:
## 1) control is a simple list, e.g. list(trace=-1)
## 2) control is a call to clmm.control
## 3) control is an empty list; list()
## all in combination with extra control arguments.

ordinal:::getCtrlArgs(clmm.control(), list(maxIter=200))
ordinal:::getCtrlArgs(list(), list(maxIter=200))
ordinal:::getCtrlArgs(list(), list(trace=-1))
ordinal:::getCtrlArgs(list(), list(trace=1))
ordinal:::getCtrlArgs(list(), list())
ordinal:::getCtrlArgs(list(maxIter=2), list())

ordinal:::getCtrlArgs(clmm.control(), list())
ordinal:::getCtrlArgs(clmm.control(maxIter=100), list(maxIter=200))
ordinal:::getCtrlArgs(clmm.control(maxIter=100), list(maxIter=200))
ordinal:::getCtrlArgs(clmm.control(), list(trace=1))
ordinal:::getCtrlArgs(clmm.control(), list(trace=-1))
ordinal:::getCtrlArgs(clmm.control(trace=1), list())
ordinal:::getCtrlArgs(clmm.control(trace=-1), list())
ordinal:::getCtrlArgs(clmm.control(trace=0), list())
## Don't specify trace twice - surprising behavior might occur: 
ordinal:::getCtrlArgs(clmm.control(trace=1), list(trace=-1))
ordinal:::getCtrlArgs(clmm.control(trace=-1), list(trace=1))
