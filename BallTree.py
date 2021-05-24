import random
import heapq

# Positive and negative infinity. In Python 3.5 can use math.inf and -math.inf
PINF =  float('inf')
NINF = -float('inf')

# node class that makes up the ballTree
class Node(object):
    
    # takes a key which is a list of points, each value
    # in the point corresponding to another dimension
    # eg: key = [(1,2,3,4), (5,3,2,1), (9,0,8,7)]
    def __init__(self, pivot = None, radius = None):
        
        self.pivot = pivot      # tuple of key and data that splits pts into
                                # left and right children
        self.radius = radius    # the dist from the pivot to the furthest pt in ball
        
        self.leftChild = None   # node with values less than pivot value
        self.rightChild = None  # node with values greater than pivot values
        
    
    def getPivotKey(self): return self.pivot[0]
    def getPivotData(self): return self.pivot[1]
    
    def __str__(self): return str(self.pivot)
    

class BallTree(object):
    
    # constructor, takes in a list of key data pairs and constructs tree of nodes
    # key's are points with numDim number of dimensions
    def __init__(self, kd, numDim):
        
        # client should not input number of dimensions <= 0
        if(numDim <= 0): raise Exception("ERROR number of dimensions must be > 0")
        
        # client must input an int as the numDim
        if not isinstance(numDim, int): raise Exception("ERROR number of dimensions must be an int")
        
        self.__root = Node()
        self.__numDim = numDim
        
        # as part of the constructor, construct the ball tree
        # starting at the root node
        self.__construct(self.__root, kd)
    
    # if pt found, return data, otherwise return None
    def find(self, ptToFind):
        
        # if the dist from the pt to the pivot of the root is greater
        # than the dist from the pivot to the furthest pt in the ball tree (radius)
        # then the pt does not exist in the ball tree so return false
        if(self.__distance(ptToFind, self.__root.getPivotKey()) > self.__root.radius):
            return None           
        
        # otherwise continue to find the pt within the tree, starting at the root
        return self.__find(ptToFind, self.__root)
    
    # invoked by find       
    def __find(self, ptToFind, curNode):
        
        # if exact match found, return data
        if(curNode.getPivotKey() == ptToFind): return curNode.getPivotData()
        
        # if both children exist and our point could be in either one,
        # recurse down both to see where it is, if it is there at all
        if(curNode.rightChild and self.__ptInBall(ptToFind, curNode.rightChild) and\
           curNode.leftChild and self.__ptInBall(ptToFind, curNode.leftChild)):
                
                leftFind = self.__find(ptToFind, curNode.leftChild)
                rightFind = self.__find(ptToFind, curNode.rightChild)
                
                # if pt found in either child, return the data
                if leftFind: return leftFind
                elif rightFind: return rightFind
                
        # if the pt could only be in the left child, recurse to find   
        elif(curNode.leftChild and self.__ptInBall(ptToFind, curNode.leftChild)):
            leftFind = self.__find(ptToFind, curNode.leftChild)
            # if found, return data
            if leftFind: return leftFind
            
        # if the pt could only be in the right child, recurse to find
        elif(curNode.rightChild and self.__ptInBall(ptToFind, curNode.rightChild)):
            # if found, return data
            rightFind = self.__find(ptToFind, curNode.rightChild)
            if rightFind: return rightFind
             
        # if we have gotten to this point, we have fully searched the tree
        # and have not found a match so return None
        return None        
    
    # invoked by find to check if a point could be within a child
    def __ptInBall(self, pt, node):
        
        # if dist <= radius: return True bc pt could possibly be found
        # in that ball
        # if dist > radius: return False bc pt couldnt possible be found
        # in that ball
        return self.__distance(node.getPivotKey(), pt) <= node.radius
    
    # returns true if leaf node, false if not
    def __isLeaf(self, cur): return not (cur.rightChild or cur.leftChild)        
    
    # input a pt and how many neighbors (N) you want to find and returns 
    # a list of N nearest neighbors to pt
    def knnFind(self, pt, N):
        
        # create a heap with N elements
        # element: (distToPt, Key, Data)
        # note: we will use heapq so we invert the values to have the implementation of a max heap
        
        # first n elements have impossible dist of NINF which
        # we will narrow down later
        # each element = tuple of (invertedDist, Key, Data)
        heap = [(NINF, None, None)]*N
        
        # heap is a list and is mutible, so mutate it
        # to include the N nearest neighbors
        # note: since heap is mutible we do not need to re-assign the var
        self.__knnFind(pt, self.__root, heap)
        
        # if there wasn't enough data to give n nearest neighbors,
        # replace 'empty' places with str 'insufficient data'
        for i in range(len(heap)):
            
            if heap[i][0] == NINF:
                heap[i] = "insufficient data"
            else:
                heap[i] = (-(heap[i][0]),) + heap[i][1:]
        
        # return updated heap
        return heap
    

    def __knnFind(self, pt, curNode, heap):

        # inverted dist
        dist = -self.__distance(curNode.getPivotKey(), pt)
        
        # worst distance is at the 0th pos of heap since
        # inverted max heap
        worstDist = heap[0][0]
        
        # if this dist from our pt to the current node is 
        # better than the worst distance in the heap, update the heap
        # dist != 0 ensures we are ignoring the point we were given
        if dist > worstDist and dist != 0:
            heapq.heappushpop(heap, (dist, curNode.getPivotKey(), curNode.getPivotData()))
        
    
        # if a right child exists and the circles formed by the
        # worstDist from our pt and the distance from our pt to the
        # pivot overlap, recurse
        if(curNode.rightChild and self.__circlesIntersect(pt, curNode.rightChild, worstDist)):

            # note: since heap is mutible we do not need to re-assign the var
            self.__knnFind(pt, curNode.rightChild, heap)
                
        # if a right child exists and the circles formed by the
        # worstDist from our pt and the distance from our pt to the
        # pivot overlap, recurse            
        if(curNode.leftChild and self.__circlesIntersect(pt, curNode.leftChild, worstDist)):
            
            # note: since heap is mutible we do not need to re-assign the var
            self.__knnFind(pt, curNode.leftChild, heap) 
        
        # if we got to this point, we have recursed to a leaf node
        # so return the heap
        return heap
        
    # circles intersect when distance btwn pivots <= sum of radii
    # invoked by knn function
    def __circlesIntersect(self, pt, cur, worstDist):
        
        distBtwnPiv = self.__distance(pt, cur.getPivotKey())
        # use neg worstDist to un-invert it 
        # note: dist was inverted to allow for max heap structure
        sumRadii = -worstDist + cur.radius                    
        
        return distBtwnPiv <= sumRadii
    
    # returns kd with no dupKeys and lets client
    # know if a value wasn't added bc of dup key
    def __noDupKeys(self, kd):
        
        keys = []
        noDup = []
        for i in range(len(kd)):
            k = kd[i][0]
            # only add kd value if k isnt dup
            if k not in keys:
                noDup+=[kd[i]]
                keys+=[k]
            else:
                print("duplicate value: ", kd, " not added")
        
        return noDup
    
    # constructs the ball tree, invoked by constructor
    def __construct(self, cur, kd):  
        
        numPts = len(kd)
        
        # if first iteration, check kd for dup keys
        # and reset it to ls without dup keys
        if cur == self.__root: 
            # if an empty kd was inserted, throw exception
            if numPts == 0: raise Exception("ERROR empty list inserted")
            kd = self.__noDupKeys(kd) 
      
        # get the dimension with the maximum spread,
        # note that dimensions start from 0
        dimGreatestSpread = self.__getDimGreatestSpread(kd)
        
        # get the pivot value
        # we will approximate the median and use that as the pivot
        numMedians = 5                               # we will perform a median of 5 to get the median 
        if(numPts < numMedians): numMedians = numPts
        medians = []                                 # to be filled with numMedian number of nodes
        
        # fill medians with numMedians random nodes in kd
        for i in range(numMedians): medians+=[kd[random.randint(0, numPts-1)]]
        
        # sort the median array on the dimension with the greatest spread
        self.__selectionSort(medians, dimGreatestSpread)
        
        # the median node is the node in the middle of the 
        # now sorted medians array, and that node will be our
        # pivot
        pivot = medians[len(medians)//2]
        
        # split based on pivot
        leftChildren = []
        rightChildren = []
        
        # split nodes based off of pivot
        # go through each pt in the list and add
        # to either left or right child cluster based off pivot
        for keyDat in kd:
            
            # data is the first element
            curKey = keyDat[0]
            pivKey = pivot[0]
            
            # if a pt doesn't have the same number of dimensions
            # as the client specified number of dimensions, throw
            # an exception
            if(len(curKey) != self.__numDim): raise Exception("ERROR point " + str(keyDat) + " has an incorrect number of dimensions")
            
            # skip over the pivot point
            if(keyDat != pivot):
            
                if(curKey[dimGreatestSpread] > pivKey[dimGreatestSpread]):
                    rightChildren += [keyDat]
                elif(curKey[dimGreatestSpread] <= pivKey[dimGreatestSpread]):
                    leftChildren += [keyDat]
        
        # at this point, leftChildren & rightChildren are lists of points (tuples) including the pt and data
        # split by the pivot value
        
        # set node attributes
        cur.pivot = pivot
        cur.radius = self.__furthestRadius(pivot, kd)
        
        # if there are left/right children pts, create the left/right
        # nodes and recurse down to construct
        
        # if leaf node end recursion by returning None
        if(leftChildren == [] and rightChildren == []): return
        
        # if there are right children and no left children
        # construct the right children
        elif(leftChildren == []):
            cur.rightChild = Node()
            self.__construct(cur.rightChild, rightChildren)
        
        # if there are left children and no right children
        # construct the left children
        elif(rightChildren == []): 
            cur.leftChild = Node()
            self.__construct(cur.leftChild, leftChildren)  
        
        # if there are both right and left children, construct both
        else:
            cur.rightChild = Node()
            cur.leftChild = Node()
            self.__construct(cur.rightChild, rightChildren)
            self.__construct(cur.leftChild, leftChildren)
             
        
        # if this was the first recursion, then set the first node
        # to be the root, the rest of the children are connected to
        # the root
        if(cur == self.__root):
            self.__root = cur
            
    def __furthestRadius(self, pivot, kd):
        
        greatestR = 0  # set the greatestRadius to lowest it could be
        
        # loop through each pt in the node
        for keyDat in kd:
            # the key is the first element
            key = keyDat[0]
            # go through each dimension
            for dim in range(self.__numDim):
                
                # get the distance btwn current point and the pivot (this is the radius)
                # the data is the first element in pivot (which is why we input pivot[0])
                dist = self.__distance(pivot[0], key)
                
                # if the distance is greater than the greatest radius
                # set the greatest radius to be this distance
                if(dist > greatestR): greatestR = dist
                    
        return greatestR   
            
            
   # return euclidean dist
   # euclidean dist = sqrt(sum from 1 to n of (qsubi - psubi)**2)
    def __distance(self, pt1, pt2):
       
        squareSums = 0.0
        
        # add together the difference squared of each dim of
        # each pt
        for dim in range(self.__numDim): squareSums += (pt1[dim] - pt2[dim])**2
        
        return squareSums**(1/2)
    
    # sorts array of nodes on value at specific dimension    
    def __selectionSort(self, kd, dim):
        
        length = len(kd)
        
        # go through selection sort
        for outer in range(length-1):
            min = outer
            for inner in range(outer+1, length):
                
                # 0 is the index of the data
                dat = kd[inner][0][dim]
                minDat = kd[min][0][dim]
                
                if dat < minDat: min = inner
            
            # swap
            kd[outer], kd[min] = kd[min], kd[outer]
    
    # return the dimension with the greatest spread
    def __getDimGreatestSpread(self,kd):
        
        # since all Nodes are assumed to have the same dim, 
        # we can just get the number of dimensions from
        # the first node in the list
        
        # create a list to store the spreads of each dim
        # so we can see which one is the greatest
        spreads = [0]*self.__numDim
        
        # find dim of greatest spread
        # for each dimension
        for dim in range(self.__numDim):
            
            # start out at most extreme values
            minVal = PINF
            maxVal = NINF
            
            # loop through each node in the inputted list
            for pt in kd:
                
                # the key is the first value in the tuple
                key = pt[0]
                
                # if that value is less than the previous minimum
                # value, set it to minVal
                if(key[dim] < minVal):
                    minVal = key[dim]
                
                # if that value is greater than the previous
                # maximum value, set it to maxVal
                if(key[dim] > maxVal):
                    maxVal = key[dim]
                    
            # set that dimension's spread to the maxVal - minVal
            spreads[dim] = maxVal-minVal
        
        # start with the max spread being the first val
        dimGreatestSpread = 0
        
        # see which dimension has the greatest spread
        for i in range(1,self.__numDim):
            if(spreads[i] > spreads[dimGreatestSpread]):
                dimGreatestSpread = i
        
        return dimGreatestSpread
    
    # print the ball tree in clusters of points           
    def pTree(self):
    
        self.__pTree(self.__root, "ROOT:  ", "")
            
    def __pTree(self, cur, kind, indent):
        
        print("\n" + indent + kind, end="")
        if cur:
            print(cur, end="")
            if cur.leftChild:
                self.__pTree(cur.leftChild, "LEFT:  ", indent + "    ")
            if cur.rightChild:
                self.__pTree(cur.rightChild, "RIGHT:  ", indent + "    ")        



