#!/bin/bash

ffmpeg -framerate 30 -i ./images/crop_res%00d.png -s 3840x2160 -c:v libx264 -vb 45000k -preset veryslow -crf 0 result.mp4 