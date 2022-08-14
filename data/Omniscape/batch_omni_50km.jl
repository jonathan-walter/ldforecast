#Batch-run omniscape using multiple .ini files

#Defaults to using 1 thread. For parallel processing, see
#Extensions -> Julia -> settings (gear icon) and scroll down to Num threads

using Omniscape #launches omniscape

#sets working directory -- REMEMBER TO CHANGE IF NEEDED
cd("/Users/kathryngrage/Box Sync/Gypsy Moth/Jeffress Spread Forecasting/JuliaScripts")

#run_omniscape("./testing/test.ini") #runs test 
#run_omniscape("./testing/test2.ini") #runs test2

run_omniscape("./Full Extent/2017_50km.ini") 
run_omniscape("./Full Extent/2018_50km.ini") 
run_omniscape("./Full Extent/2019_50km.ini") 
#run_omniscape("./Full Extent/2016_50km.ini") 
#run_omniscape("./Full Extent/2015_50km.ini") 
#run_omniscape("./Full Extent/2014_50km.ini") 
#run_omniscape("./Full Extent/2013_50km.ini") 
#run_omniscape("./Full Extent/2012_50km.ini") 
#run_omniscape("./Full Extent/2011_50km.ini") 
#run_omniscape("./Full Extent/2010_50km.ini") 
#run_omniscape("./Full Extent/2009_50km.ini") 
#run_omniscape("./Full Extent/2008_50km.ini") 
#run_omniscape("./Full Extent/2007_50km.ini") 
#run_omniscape("./Full Extent/2006_50km.ini") 
#run_omniscape("./Full Extent/2005_50km.ini") 
#run_omniscape("./Full Extent/2004_50km.ini") 
#run_omniscape("./Full Extent/2003_50km.ini") 
#run_omniscape("./Full Extent/2002_50km.ini") 
#run_omniscape("./Full Extent/2001_50km.ini") 