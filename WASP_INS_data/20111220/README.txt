The following is a description of the experiments in the INS data files in this folder.

Using NED coordinates, while strapped to the trolley, 315 is to the left (west), 220 is to the right (east),
and they are separated by 1.145m.

Most of these experiments (except maybe the first two) are in the .mov file, which is filmed looking at the origin from the 'right'


Exp	Type		w220 file	w315 file	WASP file	Description

1	Trolley		1		1		1		0m to 24m straight line and stop, on RHS of centre line, trolley front wheels start and stop on the line		
2	Trolley		2		2		2		24m to 0m straight line and stop, on LHS of centre line, trolley starts back wheels on the line and stops front wheels on the line
3	Trolley		3		3		3		At 2m point, starting facing to the origin, 1 full rotation right then 1 full rotation left, stopping at each 90 degrees
4	Trolley		4		4		4		0m to 50m (passing outside of convex hull), RHS of centre line, trolley front wheels start and stop on line, rotate 90 degrees and return to start along same tracks
5	Trolley		5		5		5		Around the outside of concrete circle, starting at crack going anti-clockwise one revolution, then going backwards clockwise to finish at the starting position. Possible magnetic interference, especially around halfway of each revolution from two railroad tracks
6	Trolley		6		6		6		Boxing the L shaped path, anticlockwise, see map, possible magnetic interference at about half way, which is also outside of the convex hull, ending where it started
7	Trolley		7		7		NA		Calibrating the sensors for orientation, at origin facing towards 50m, small adjustment made mid way (use 2nd half of log), then a full revolution to try and calibrate magnetic vector
8	Sathyan's trial	NA		NA		7		NA
9	Backpack	NA		8		8		Repeating Exp 6 around the path again at a brisk walk, but with a single node in a backpack, stopping abruptly at 90 degrees turns, again ending where it started
10	Backpack	NA		9		9		Repeating Exp 5 around the circle with the backpack, brisk walk, but this time orientation is changed before retracing the circle (it's not walking backwards)
11	Backpack	NA		10		10		Running Exp 6 with backpack, again stopping at 90 degrees corners
12	Backpack	NA		11		11		Running around randomly, on the grass and around the circle, including shuffling sideways on the grass and dropping to one knee
13	Static		NA		12		12		Inside the convex hull
14	Static		NA		13		13		Inside	
15	Static		NA		14		14		Inside
16	Static		NA		15		15		Outside the convex hull


WASP Data Format:
1. 	8 bit overflowing sequence number
2. 	Node ID (220)
3. 	Valid bit
4. 	Timestamp
5. 	x
6. 	y
7. 	z
8. 	error (std residuals?)
9. 	Node ID (315)
10. 	Valid bit
11. 	Timestamp
12. 	x
13. 	y
14. 	z
15. 	error (std residuals?)


Survey Data Format:
Ignore the first line, then headers are:
1.	Node	
2.	x
3.	y
4.	z	
5.	valid


filemapping.csv
1.	wsd (IMU) data file
2. 	associated WASP data file
3. 	which col of data to use out of the WASP data file
4-6.	vector towards the 'front' of the WASP unit on that run (for plotting only)
7-9.	calibrated magnetometer bias for that node


