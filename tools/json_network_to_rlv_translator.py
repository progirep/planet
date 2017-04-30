#!/usr/bin/env python2
import os, sys
import json, operator
from pprint import pprint

# Read input file
if len(sys.argv)<2:
    print >>sys.stderr, "Error: Expected JSON input file name!"
    
with open(sys.argv[1]) as data_file:    
    data = json.load(data_file)
data = data["layer"]

# Neuron lookup table: layer--->Neuron names
neurons = {}


# Function for recursively processing layers
def recurseProcessLayer(dataLineName,dataLineWidth):
    '''Returns the neuron names. DataLineWidth may be unknown.'''

    # DAG Caching
    if dataLineName in neurons:
        return neurons[dataLineName]
        
    # Search for the producer of this layer
    for layer in data:
        if dataLineName in layer["top"]:

            # This is the layer to be processed
            # ---> Proceed according to type.
            if layer["type"]=="Split":
                assert len(layer["bottom"])==1
                return recurseProcessLayer(layer["bottom"][0],dataLineWidth)
            elif layer["type"]=="InnerProduct":
                
                # Get output dimension
                outputLineWidth = layer["inner_product_param"]["num_output"]
                if dataLineWidth!=None:
                    assert outputLineWidth==dataLineWidth
                
                # Now get input dimension and blob data
                inputLineWidth = None
                biasWeights = None
                inputWeights = None
                for blob in layer["blobs"]:
                    if len(blob["shape"]["dim"])==1:
                        biasWeights = blob["data"]
                    elif len(blob["shape"]["dim"])==2:
                        inputLineWidth = blob["shape"]["dim"][1]
                        inputWeights = blob["data"]
                        assert blob["shape"]["dim"][0]==outputLineWidth
                    else:
                        assert False
                assert inputLineWidth != None
                
                # Now get input
                assert len(layer["bottom"])==1
                inputNeurons = recurseProcessLayer(layer["bottom"][0],inputLineWidth)
                
                # Produce outputs
                outputNeurons = []
                for i in range(0,outputLineWidth):
                    outputNeurons.append(dataLineName+"X"+str(i))
                    sys.stdout.write("Linear "+dataLineName+"X"+str(i)+" "+str(biasWeights[i]))
                    for j in range(0,len(inputNeurons)):
                        sys.stdout.write(" "+str(inputWeights[i*inputLineWidth+j])+" "+str(inputNeurons[j]))
                    sys.stdout.write("\n")
                
                neurons[dataLineName] = outputNeurons
                return outputNeurons
                
                
            elif layer["type"]=="ReLU":
            
                # RELU: Get Input
                assert len(layer["bottom"])==1
                inputNeurons = recurseProcessLayer(layer["bottom"][0],dataLineWidth)
                
                # Produce outputs
                outputNeurons = []
                for i in range(0,len(inputNeurons)):
                    outputNeurons.append(dataLineName+"X"+str(i))
                    sys.stdout.write("ReLU "+dataLineName+"X"+str(i)+" 0.0 1.0 "+str(inputNeurons[i])+"\n")
                
                neurons[dataLineName] = outputNeurons
                return outputNeurons
                
                
                
            elif layer["type"]=="Convolution":
                
                # Convolution
                
                # ---> Load weights and biasses
                dimWeights = None
                weights = None
                biasses = None
                for blob in layer["blobs"]:
                    size = blob["shape"]["dim"]
                    params = blob["data"]
                    if len(size)==1:
                        biasses = params
                    else:
                        dimWeights = size
                        weights = params
                        
                # ---> Other parameters
                #--------> Read Stride
                if "stride" in layer["convolution_param"]:
                    stride = layer["convolution_param"]["stride"]
                    assert len(stride)==1
                    stride = stride + stride
                else:
                    if "stride_w" in layer["convolution_param"]:
                        if "stride_h" in layer["convolution_param"]:
                            stride = [layer["convolution_param"]["stride_w"],layer["convolution_param"]["stride_h"]]
                        else:
                            stride = [layer["convolution_param"]["stride_w"],1]
                    else:
                        if "stride_h" in layer["convolution_param"]:
                            stride = [1,layer["convolution_param"]["stride_h"]]
                        else:
                            stride = [1,1]
                
                #--------> Read Kernel Size
                if "kernel_size" in layer["convolution_param"]:
                    kernel_size = layer["convolution_param"]["kernel_size"]
                    assert len(kernel_size)==1
                    kernel_size = kernel_size + kernel_size
                else:
                    kernel_size = [layer["convolution_param"]["kernel_w"],layer["convolution_param"]["kernel_h"]]

                #--------> Read PAD
                if "pad" in layer["convolution_param"]:
                    padding = layer["convolution_param"]["pad"]
                    assert len(padding)==1
                    padding = padding + padding
                else:
                    if "pad_w" in layer["convolution_param"]:
                        if "pad_h" in layer["convolution_param"]:
                            padding = [layer["convolution_param"]["pad_w"],layer["convolution_param"]["pad_h"]]
                        else:
                            padding = [layer["convolution_param"]["pad_w"],0]
                    else:
                        if "pad_h" in layer["convolution_param"]:
                            padding = [0,layer["convolution_param"]["pad_h"]]
                        else:
                            padding = [0,0]

                num_output = layer["convolution_param"]["num_output"]
                num_input_channels = dimWeights[1]
                
                # Rest is unimplemented for the time being.
                # assert dimWeights[0]==1 #---> This is for the *outgoing* num_outputs
                # assert dimWeights[1]==1 #----> This is for the incoming colors or features
                
                # Check for some unsupported features
                if "bias_term" in layer["convolution_param"]:
                    print >>sys.stderr, "Error: Only the default 'bias_term' value is supported for convolution layers"
                    sys.exit(1)

                if "group" in layer["convolution_param"]:
                    print >>sys.stderr, "Error: Only the default 'group' value is supported for convolution layers"
                    sys.exit(1)
                
                # ---> Read input
                inputNeurons = recurseProcessLayer(layer["bottom"][0],None)
                
                # ---> Unflatten weights
                def unflatten(neurons,remainingDimensions):
                    if len(remainingDimensions)==1:
                        return (neurons[0:remainingDimensions[0]],neurons[remainingDimensions[0]:])
                    else:
                        res = []
                        for a in range(0,remainingDimensions[0]):
                            (d,neurons) = unflatten(neurons,remainingDimensions[1:])
                            res.append(d)
                        return (res,neurons)
                
                (unflattenedWeights,rest) = unflatten(weights,dimWeights)
                assert rest==[]
                
                # Compute convolution
                resultingNeurons = []
                for i in range(0,num_output):
                    ysize = len(inputNeurons[0])
                    xsize = len(inputNeurons[0][0])
                    thisBlock = []
                    for y in xrange(-1*padding[1],ysize-kernel_size[1]+1+padding[1],stride[1]):
                        thisLine = []                    
                        for x in xrange(-1*padding[0],xsize-kernel_size[0]+1+padding[0],stride[0]):
                            thisLine.append(dataLineName+"X"+str(i)+"X"+str(x)+"X"+str(y))
                            localInputs = []
                            for c in range(0,num_input_channels):
                                for b in range(0,kernel_size[1]):
                                    for a in range(0,kernel_size[0]):
                                        if y+b>=0 and y+b<len(inputNeurons[c]):
                                            if x+a>=0 and x+a<len(inputNeurons[c][y+b]):
                                                localInputs.append(str(unflattenedWeights[i][c][b][a])+" "+inputNeurons[c][y+b][x+a]) 
                            sys.stdout.write("Linear "+thisLine[-1]+" "+str(biasses[i])+" "+" ".join(localInputs)+"\n")
                        thisBlock.append(thisLine)
                    resultingNeurons.append(thisBlock)
                    
                return resultingNeurons

            elif layer["type"]=="HDF5Data":
            
                # Input layer
                assert dataLineWidth!=None
                
                outputNeurons = []
                for i in range(0,dataLineWidth):
                    outputNeurons.append("inX"+str(i))
                    sys.stdout.write("Input inX"+str(i)+"\n")

                neurons[dataLineName] = outputNeurons
                return outputNeurons
                
            elif layer["type"]=="Data":
            
                # Input layer
                assert dataLineWidth!=None
                
                # We assume normalization by a factor of 1/256.0 here.
                assert layer['transform_param']['scale'] == 0.00390625

                outputNeurons = []
                for i in range(0,dataLineWidth):
                    outputNeurons.append("inX"+str(i))
                    sys.stdout.write("Input inX"+str(i)+"\n")

                neurons[dataLineName] = outputNeurons
                return outputNeurons
                
            elif layer["type"]=="Reshape":
            
                # Reshape layer
                outputDimension = layer['reshape_param']['shape']['dim']
                # print layer
                assert outputDimension[0]==-1 # The first dimension is always the sample points
                
                nofInputs = 1
                for a in outputDimension[1:]:
                    nofInputs *= a
                inputNeurons = recurseProcessLayer(layer["bottom"][0],nofInputs)

                
                # Ok, first flatten input Neurons
                def flat(neurons):
                    if type(neurons)==str or type(neurons)==unicode:
                        return [neurons]
                    else:
                        l = []
                        for a in neurons:
                            l.extend(flat(a))
                        return l
                
                # Ok, first flatten input Neurons
                #def flat(neurons):
                #    if type(neurons)==str or type(neurons)==unicode or type(neurons)==int:
                #        return [neurons]
                #    else:
                #        l = []
                #        for a in neurons:
                #            l.append(flat(a))
                #        return [a for j in zip(*l) for a in j ]

                flattenedNeurons = flat(inputNeurons)
                
                # Ok, now reshape
                def unflatten(neurons,remainingDimensions):
                    if len(remainingDimensions)==1:
                        return (neurons[0:remainingDimensions[0]],neurons[remainingDimensions[0]:])
                    else:
                        res = []
                        for a in range(0,remainingDimensions[0]):
                            (d,neurons) = unflatten(neurons,remainingDimensions[1:])
                            res.append(d)
                        return (res,neurons)
                
                # Ok, now reshape
                #def unflatten(neurons,selection,dimensions):
                #    if len(selection)==len(dimensions):
                #        index = 0
                #        factor = 1
                #        for i in range(0,len(selection)):
                #            index += factor*selection[i]
                #            factor *= dimensions[i]
                #        return neurons[index]
                #    else:
                #        res = []
                #        for a in range(0,dimensions[len(selection)]):
                #            res.append(unflatten(neurons,selection+[a],dimensions))
                #        return res

                
                (unflattenedNeurons,rest) = unflatten(flattenedNeurons,outputDimension[1:])
                assert rest==[]
                # print outputDimension
                # print inputNeurons
                # print flattenedNeurons
                # print "-----UF:->",unflattenedNeurons
                assert len(flattenedNeurons)==reduce(operator.mul, outputDimension[1:], 1)
                              
                return unflattenedNeurons

            elif layer["type"]=="Pooling":
            
                # Reshape layer
                inputNeurons = recurseProcessLayer(layer["bottom"][0],None)
                
                assert layer["pooling_param"]["pool"]==0 # Must be a MAXPOOL (for the time being)
                
                # ---> Other parameters
                #--------> Read Stride
                if "stride" in layer["pooling_param"]:
                    # Why is the "stride" here an int, but for the Convolution layer is a list?
                    stride = layer["pooling_param"]["stride"]
                    if isinstance(stride,list):
                        assert len(stride)==1
                        stride = stride + stride
                    else:
                        stride = [stride,stride]
                else:
                    if "stride_w" in layer["pooling_param"]:
                        if "stride_h" in layer["pooling_param"]:
                            stride = [layer["pooling_param"]["stride_w"],layer["pooling_param"]["stride_h"]]
                        else:
                            stride = [layer["pooling_param"]["stride_w"],1]
                    else:
                        if "stride_h" in layer["pooling_param"]:
                            stride = [1,layer["pooling_param"]["stride_h"]]
                        else:
                            stride = [1,1]
                
                #--------> Read Kernel Size
                if "kernel_size" in layer["pooling_param"]:
                    # Why is the "kernel_size" here an int, but for the Convolution layer is a list?
                    kernel_size = layer["pooling_param"]["kernel_size"]
                    if isinstance(kernel_size,list):
                        assert len(kernel_size)==1
                        kernel_size = kernel_size + kernel_size
                    else:
                        kernel_size = [kernel_size,kernel_size]
                else:
                    kernel_size = [layer["pooling_param"]["kernel_w"],layer["pooling_param"]["kernel_h"]]

                #--------> Read PAD
                if "pad" in layer["pooling_param"]:
                    padding = layer["pooling_param"]["pad"]
                    assert len(padding)==1
                    padding = padding + padding
                else:
                    if "pad_w" in layer["pooling_param"]:
                        if "pad_h" in layer["pooling_param"]:
                            padding = [layer["pooling_param"]["pad_w"],layer["pooling_param"]["pad_h"]]
                        else:
                            padding = [layer["pooling_param"]["pad_w"],0]
                    else:
                        if "pad_h" in layer["pooling_param"]:
                            padding = [0,layer["pooling_param"]["pad_h"]]
                        else:
                            padding = [0,0]

                
                # Here, we assume that the "inputNeurons" array is three-dimensional:
                # - color channel
                # - X channel
                # - Y channel
                resultingNeurons = []
                for i,channel in enumerate(inputNeurons):
                    ysize = len(channel)
                    xsize = len(channel[0])
                    thisBlock = []
                    for y in xrange(-1*padding[1],ysize-kernel_size[1]+1+padding[1],stride[1]):
                        thisLine = []                    
                        for x in xrange(-1*padding[0],xsize-kernel_size[0]+1+padding[0],stride[0]):
                            thisLine.append(dataLineName+"X"+str(i)+"X"+str(x)+"X"+str(y))
                            localInputs = []
                            for b in range(0,kernel_size[1]):
                                for a in range(0,kernel_size[0]):
                                    if y+b>=0 and y+b<len(channel):
                                        if x+a>=0 and x+a<len(channel[y+b]):
                                            localInputs.append(channel[y+b][x+a]) 
                            sys.stdout.write("MaxPool "+thisLine[-1]+" "+" ".join(localInputs)+"\n")
                        thisBlock.append(thisLine)
                    resultingNeurons.append(thisBlock)
                return resultingNeurons                    
                    
            else:
                raise RuntimeError("Unsupported Layer Type: "+layer["type"])
            
    raise "Error: Data Line "+dataLineName+" not found."


# Process the Accuracy layer
foundAccurracyLayer = False
for layer in data:
    if layer['type'] == 'Accuracy':
        foundAccurracyLayer = True
        outputs = recurseProcessLayer(layer["bottom"][0],None)
        for i in range(0,len(outputs)):
            sys.stdout.write("Linear outX"+str(i)+" 0.0 1.0 "+outputs[i]+"\n")
                
if not foundAccurracyLayer:
    print >>sys.stderr, "Warning: No 'Accuracy' layer found, hence nothing was translated."
