system("make");
system("./gquad.x > gquad.txt");
system("./gquad2.x > gquad2.txt");
system("./trappx.x > trappx.txt");
system("./trappx2.x > trappx2.txt");

system("nohup matlab -nosplash -nodisplay < appx_error.m > output.txt");
system("mv ErrorPlot.png ../Report");
