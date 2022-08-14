#Batch-run omniscape using multiple .ini files

#Defaults to using 1 thread. For parallel processing, see
#Extensions -> Julia -> settings (gear icon) and scroll down to Num threads

using Omniscape #launches omniscape

#sets working directory -- REMEMBER TO CHANGE IF NEEDED
cd("/Users/kathryngrage/Box Sync/Gypsy Moth/Jeffress Spread Forecasting/JuliaScripts")

#run_omniscape("./testing/test.ini") #runs test 
#run_omniscape("./testing/test2.ini") #runs test2

run_omniscape("./Full Extent/2017_10km.ini") 
run_omniscape("./Full Extent/2018_10km.ini") 
run_omniscape("./Full Extent/2019_10km.ini") 
#run_omniscape("./Full Extent/2016_10km.ini") 
#run_omniscape("./Full Extent/2015_10km.ini") 
#run_omniscape("./Full Extent/2014_10km.ini") 
#run_omniscape("./Full Extent/2013_10km.ini") 
#run_omniscape("./Full Extent/2012_10km.ini") 
#run_omniscape("./Full Extent/2011_10km.ini") 
#run_omniscape("./Full Extent/2010_10km.ini") 
#run_omniscape("./Full Extent/2009_10km.ini") 
#run_omniscape("./Full Extent/2008_10km.ini") 
#run_omniscape("./Full Extent/2007_10km.ini") 
#run_omniscape("./Full Extent/2006_10km.ini") 
#run_omniscape("./Full Extent/2005_10km.ini") 
#run_omniscape("./Full Extent/2004_10km.ini") 
#run_omniscape("./Full Extent/2003_10km.ini") 
#run_omniscape("./Full Extent/2002_10km.ini") 
#run_omniscape("./Full Extent/2001_10km.ini") 