@echo off

ffmpeg -r 2 -start_number 0 -i "sim\iter %%06d.jpg" -c:v libx264 -crf 24 -r 24 -pix_fmt yuv420p out.mp4