% define the Room dimensions
roomDimensions = [4 4 2.5];
% define the coordinates of the receiver and transmitter
receiverCoord = [2 1 1.8];
sourceCoord = [3 1 1.8];
% plot the room, receiver, and transmitter
h = figure;
plotRoom(roomDimensions,receiverCoord,sourceCoord,h)

