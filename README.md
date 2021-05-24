# **BallTree Class**
BallTree class takes in a list of key data pairs and an integer for the number of dimensions you will be using

## Introduction
A Ball Tree is a data structure useful for working in multi-dimensional space. It is made up of nodes which have a key (point), data, a left child, and a right child. This structure is particularly useful for nearest neighbor search in many dimensions. 

## This Implementation
My Ball Tree's left and right children are single nodes though there are other implmentations where left and right children are groups of nodes. 

My class supports find and knnFind.

## Construction
In this implementation, construction of the Ball Tree takes place when the tree is instantiated. The object is passed a list of tuples as well as the number of dimensions you will be working in. The tuple includes the key and data pair. 

The tree is constructed by going through each key data pair, choosing a pivot point, splitting the points on the pivot, and recursively constructing each side until all points are in the tree. 

## Implementation Examples
### Construction/Initializing Tree

`
kd = [((1,2,3), "A"), ((4,5,6), "B"), ((1,4,5), "C"), ((4,2,2), "D"), ((1,1,5), "E")]
numDimensions = 3
t = BallTree.BallTree(kd, numDimensions)
`

### Printing the Tree
`t.pTree()`
`ROOT:  ((4, 5, 6), 'B')
    LEFT:  ((1, 1, 5), 'E')
       LEFT:  ((1, 4, 5), 'C')
            LEFT:  ((1, 2, 3), 'A')
        RIGHT:  ((4, 2, 2), 'D')`

### Find: 
`pt = (4,5,6)  
t.find(pt)  `

`B`

### KnnFind: 
`
pt = (1,2,3)  
N = 2  
t.knnFind(pt, N)  `

`[(2.8284271247461903, (1, 4, 5), 'C'), (2.23606797749979, (1, 1, 5), 'E')]`
