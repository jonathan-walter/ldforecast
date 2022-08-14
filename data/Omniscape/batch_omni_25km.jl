#Batch-run omniscape using multiple .ini files

#Defaults to using 1 thread. For parallel processing, see
#Extensions -> Julia -> settings (gear icon) and scroll down to Num threads

using Omniscape #launches omniscape

#sets working directory -- REMEMBER TO CHANGE IF NEEDED
cd("/Users/kathryngrage/Box Sync/Gypsy Moth/Jeffress Spread Forecasting/JuliaScripts")

#run_omniscape("./testing/test.ini") #runs test 
#run_omniscape("./testing/test2.ini") #runs test2

run_omniscape("./Full Extent/2017_25km.ini") 
run_omniscape("./Full Extent/2018_25km.ini") 
run_omniscape("./Full Extent/2019_25km.ini") 
#run_omniscape("./Full Extent/2016_25km.ini") 
#run_omniscape("./Full Extent/2015_25km.ini") 
#run_omniscape("./Full Extent/2014_25km.ini") 
#run_omniscape("./Full Extent/2013_25km.ini") 
#run_omniscape("./Full Extent/2012_25km.ini") 
#run_omniscape("./Full Extent/2011_25km.ini") 
#run_omniscape("./Full Extent/2010_25km.ini") 
#run_omniscape("./Full Extent/2009_25km.ini") 
#run_omniscape("./Full Extent/2008_25km.ini") 
#run_omniscape("./Full Extent/2007_25km.ini") 
#run_omniscape("./Full Extent/2006_25km.ini") 
#run_omniscape("./Full Extent/2005_25km.ini") 
#run_omniscape("./Full Extent/2004_25km.ini") 
#run_omniscape("./Full Extent/2003_25km.ini") 
#run_omniscape("./Full Extent/2002_25km.ini") 
#run_omniscape("./Full Extent/2001_25km.ini") 