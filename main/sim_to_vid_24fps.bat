@echo off

ffmpeg -y -r 24 -start_number 0 -i "sim\waves %%06d.jpg" -vf scale=-1:1080:flags=neighbor -c:v libx264 -crf 24 -r 24 out.mp4