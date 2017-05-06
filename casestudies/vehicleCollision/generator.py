#!/usr/bin/env python2
#
# Testing pyBrain
#
import os, pygame, pygame.locals, sys, random, math

# CONSTANTS
XSIZE = 480
YSIZE = 480

# OBSTACLE/ROBOT RADIUS
OBSTACLE_INNER_RADIUS = 10
ROBOT_INNER_RADIUS = 10
JOINT_OUTER_RADIUS = 30


# Prepare screen
screen = pygame.display.set_mode((XSIZE,YSIZE))
pygame.display.set_caption('Collision Yes/No Visualizer')
clock = pygame.time.Clock()
isPaused = False
screenBuffer = pygame.Surface(screen.get_size())
screenBuffer = screenBuffer.convert()

# Buffer for results
COLLISION_CASES = []
NON_COLLISION_CASES = []

# distance computation
def computeDistance(a,b):
    diffA = a[0]-b[0]
    diffB = a[1]-b[1]
    return math.sqrt(diffA*diffA+diffB*diffB)


outFile = open("collisions.csv","a")

# Main loop
nofExamplesWritten = 0
while 1:

    currentPos = (XSIZE/2,YSIZE*5/6)
    currentSpeed = 0.2
    currentDirection = math.pi
    directionChange = random.random()*0.0012-0.0006
    
    obstaclePos = (random.random()*XSIZE,random.random()*YSIZE)
    obstacleSpeed = 0.1+random.random()*0.2
    obstacleDirection = (0.4+0.2*random.random())*math.pi
    obstacleDirectionChange = random.random()*0.0012-0.0006
    
    
    case = (obstaclePos[0]/XSIZE,obstaclePos[1]/YSIZE,obstacleSpeed,obstacleDirection/math.pi,obstacleDirectionChange/0.0020,directionChange/0.0020)
    
    
    collision = 0    # 1: Maybe, 2: Sure
    
    screenBuffer.fill((64, 64, 64))
    for i in range(0,3000):
        if (i % 10) == 0:
            if pygame.display.get_active():
                if collision==0:
                    collisionColor=0
                elif collision==1:
                    collisionColor=128
                elif collision==2:
                    collisionColor=255
                pygame.draw.circle(screenBuffer,(255,255,collisionColor),(int(currentPos[0]),int(currentPos[1])),ROBOT_INNER_RADIUS)
                pygame.draw.circle(screenBuffer,(255,collisionColor,255),(int(obstaclePos[0]),int(obstaclePos[1])),OBSTACLE_INNER_RADIUS)
        currentPos = (math.sin(currentDirection)*currentSpeed+currentPos[0],math.cos(currentDirection)*currentSpeed+currentPos[1])
        currentDirection += directionChange
        obstaclePos = (math.sin(obstacleDirection)*obstacleSpeed+obstaclePos[0],math.cos(obstacleDirection)*obstacleSpeed+obstaclePos[1])
        obstacleDirection += obstacleDirectionChange
        
        # Distance
        distance = computeDistance(currentPos,obstaclePos)
        if distance<JOINT_OUTER_RADIUS:
            collision = max(collision,1)
        if distance<OBSTACLE_INNER_RADIUS+ROBOT_INNER_RADIUS:
            collision = 2
        

    if pygame.display.get_active():
        clock.tick(1)
    screen.blit(screenBuffer,(0,0))
    pygame.display.flip()
    

    # Sort into the collision cases / non-collision cases    
    if collision==0:
        NON_COLLISION_CASES.append(case)
    elif collision==2:
        COLLISION_CASES.append(case)
        
        
    # Lists are too long? Then Truncate
    if len(NON_COLLISION_CASES)>100:
        NON_COLLISION_CASES = NON_COLLISION_CASES[0:100]
    if len(COLLISION_CASES)>100:
        COLLISION_CASES = COLLISION_CASES[0:100]
        
    # We have both collision and non-collision cases? Then flush them out!
    if len(NON_COLLISION_CASES)>0 and len(COLLISION_CASES)>0:
        outFile.write(",".join([str(a) for a in COLLISION_CASES[0]])+",1\n")
        outFile.write(",".join([str(a) for a in NON_COLLISION_CASES[0]])+",0\n")
        outFile.flush()
        nofExamplesWritten += 2
        NON_COLLISION_CASES = NON_COLLISION_CASES[1:]
        COLLISION_CASES = COLLISION_CASES[1:]
        
        if (nofExamplesWritten % 100)==0:
            sys.stdout.write(str(nofExamplesWritten)+"...")
            sys.stdout.flush()

    
    #  Event loop
    for event in pygame.event.get():
        
        if event.type == pygame.locals.QUIT or (event.type == pygame.locals.KEYDOWN and event.key == pygame.locals.K_ESCAPE):
            # from pybrain.tools.xml.networkwriter import NetworkWriter
            # NetworkWriter.writeToFile(net, 'finalNet.xml')            
            sys.exit(0)
        if (event.type == pygame.locals.KEYDOWN and event.key == pygame.locals.K_SPACE):
            isPaused = not isPaused


