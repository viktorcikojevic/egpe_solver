# !/bin/bash
rm out.mp4
rm out.gif
ffmpeg -framerate 25 -start_number 0 -i  2D_snapshot_%d_xy.png -c:v libx264 -r 40  out.mp4
ffmpeg -i out.mp4 out.gif
gifsicle --resize 512x512 --colors 8 out.gif > smaller.gif
mv smaller.gif out.gif
