import os

for i in range(0, 1500):
    print(i)
    os.system("convert res" + str(i) + ".png -crop 5262x2960+2041+1080 crop_res" + str(i) + ".png")
