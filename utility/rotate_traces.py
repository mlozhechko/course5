import os

# config
time = 50
fps = 30
min_angle = 0
max_angle = 1
initial_angle = -0.03
angle_around_x = 0.4

delta_angle = max_angle - min_angle
step_angle = delta_angle / (time * fps)

i = 0
angle_i = min_angle
while i < fps * time:
    print(str(i) + " step of " + str(fps * time))
    os.system("../course -x 2400 -y 1800 -j 16 -f ../test2.vtk -d result.vti -X " + str(angle_around_x) + " -I " + str(initial_angle) + " -Y " + str(angle_i))
    os.system("pvpython screen.py " + str(i))
    i += 1
    angle_i += step_angle
