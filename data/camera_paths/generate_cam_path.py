import math
import csv
import sys

def gen_straight_path2(filename, initialCam, deltaMove, deltaLook): 
    cx = initialCam[0]
    cy = initialCam[1]
    cz = initialCam[2]

    lookx = initialCam[3]
    looky = initialCam[4]
    lookz = initialCam[5]

    navigable_positions = []
    for i in range(450):
        spot = []

        spot.append(cx)
        spot.append(cy)
        spot.append(cz)
        spot.append(lookx)
        spot.append(looky)
        spot.append(lookz)
                                                           
        navigable_positions.append(spot)

        cx -= deltaMove
        lookz += deltaLook

    with open(filename,'w') as f:
        f.writelines(' '.join(str(j) for j in i) +'\n' for i in navigable_positions)

def main():
    scene = sys.argv[1]
    scenePos = open(scene + ".txt", "r")
    data = [[float(i) for i in line.split()] for line in scenePos]

    initCam = [data[0][0], data[0][1], data[0][2], data[0][4], data[0][5], data[0][6]]
    deltaMove = 1.0/100.0
    deltaLook = 1.0/250.0
    gen_straight_path2("room1-camera2.txt", initCam, deltaMove, deltaLook)

    #initCam = [0.5, 0.5, 0.5, -3, 0.5, -0.4]
    #deltaMove = 1.0/100.0
    #deltaLook = 1.0/250.0
    #gen_straight_path2("room1-camera1.txt", initCam, deltaMove, deltaLook)

if __name__ == '__main__':
    main()
