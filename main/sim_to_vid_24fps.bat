@echo off

ffmpeg -y -r 24 -start_number 0 -i "sim\iter %%06d.png" -vf scale=-1:1080:flags=neighbor -c:v libx264 -crf 32 -r 24 out.mp4