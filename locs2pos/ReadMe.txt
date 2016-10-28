Usage:
locs2txt convert the cluster location file in "locs" format into plain text.
The command
	locs2txt fileA fileB
will read the locs file fileA, convert it into _pos.txt format, and output it as fileB.

locs2txt_all is a shell script, which convert all the "locs" files in a directory into 
plain text by calling the above tool locs2txt.
The command
	./locs2txt_all.sh dirA dirB
will onvert all the locs files in the directory dirA into _pos.txt and output them in dirB.