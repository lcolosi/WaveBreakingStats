# Code Development Log 
Here, I will document the changes made to the wave breaking statistic scripts. Additionally, questions/comments concerning data processing, potential mistakes, and possible improvements I can make to the code will be included.  

----
### Main.m 
1. Moved [documentation](https://github.com/lcolosi/WaveBreakingStats/tree/main/src) to the Github repository.
2. Consolidated plug-and-play parameters to the beginning of the script. 
3. 

### TrackSteady.m 
1. Removed previous approach for averaging heading (adding 360 to values close to wrapping point). Now, the mean angle is computed for the average cosine and sine values of the heading angle. The standard deviation of heading is computed using the Yamartino method. The mean and standard deviation of roll and pitch are computed using simply the mean and standard deviation function built-into matlab (the roll and pitch don't approach the wrapping point).

2. Changed code for detecting abrupt changes in angles. Recall that abrupt changes wereoriginaly deteched using two approaches in the code: 

> The jth roll/pitch/heading observation deviates from the its mean by more than t times its standard deviation plus 1 (for roll and pitch) or 2 (for heading).

> The jth roll/pitch/heading observation 7 seconds away (when the jth observation is within 7 seconds, use the time index at the beginning/end of the track) deviates from its mean by more than t times its standard deviation plus 1 (for roll and pitch) or(for heading).

There seemed to be two mistakes here: 
    1. There is no reason I can tell for adding 1 or 2 to `t*\sigma`

**Questions**
1. What does the time value represent in the EO files? Is it seconds of the week? 

**Comment** 
1. Instead of using the standard deviation of the roll, pitch, and heading to figure out when the plane is turning off the flight line, the gradient of the roll pitch and heading (computed over 5-10 second intervals) may be used.  

