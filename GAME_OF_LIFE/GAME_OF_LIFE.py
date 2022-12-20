#!/usr/bin/env python
# coding: utf-8

def board_print(board,size):#define board printing function
    for x in range(size+2):
        print("-",end="")#print edges
    print("\n",end="")
    for x in range(size):
        print("|",end="")
        for y in range(size):
            if board[x][y] == 1:#if position in board == 1
                print("X",end="")#print X
            else:
                print(" ",end="")#print " " if position in board == 0
        print("|\n",end="")
    for x in range(size+2):
        print("-",end="")
    print("\n",end="")
    
def compute_neighbours(board,size,row,column):#count neighbouring cells
    x=row-1#assign neighbouring positions (rows and columns) 
    x1=row
    x2=row+1
    y=column-1
    y1=column
    y2=column+1
    rows=[x,x1,x2]
    columns=[y,y1,y2]
    positions=[]
    neighbour_count=0
    for x in rows:
        for y in columns:
            positions.append([x,y])#create a list of 2-element lists (position of neighbours)
            for i in positions:
                if i[0]<0 or i[1]<0 or i[0]>size-1 or i[1]>size-1:#if any row or column is < 0 
                    positions.remove(i)#delete it from the list
    positions.remove([x1,y1])#remove cell whose neighbours we are looking for
    for i in positions:
        neighbour_count+=board[i[0]][i[1]]#count how many living (+1) neighbours cell has (dead neighbours add 0)
    return neighbour_count

def next_generation(board,size):
    next_board = [[0 for _ in range(size)] for _ in range(size)]
    for x in range(size):
        for y in range(size):
            neighbours = compute_neighbours(board,size,x,y)#compute neighbours for a given cell
            if board[x][y] == 1 and (neighbours < 2 or neighbours > 3):#if neighbours <2 or >3
                next_board[x][y] = 0 #set value to 0 in this position for the next matrix (cell dies)
            elif board[x][y] == 1 and (neighbours==2 or neighbours==3):
                next_board[x][y] = 1 #cells that keep living
            elif board[x][y] == 0 and neighbours == 3:
                next_board[x][y] = 1 #cells that revive
            else:
                next_board[x][y]=board[x][y]#other situations remain the same
    return next_board

#import packages
import random
from time import sleep
size=10#define size
board = [[random.randint(0,1) for _ in range(size)] for _ in range(size)]#create starting random board
print("Next generation every 2 seconds.\nPress Ctrl+C to end simulation.")
while True:
    board_print(board,size)#print first board
    sleep(2)#wait 2 seconds
    board= \
    next_generation(board,size)#print next generation board







