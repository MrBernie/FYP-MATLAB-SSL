function plotRoom(roomDimensions,receiverCoord,sourceCoord,figHandle)
% PLOTROOM Plot room, transmitter and receiver

figure(figHandle)
X = [0;roomDimensions(1);roomDimensions(1);0;0];
Y = [0;0;roomDimensions(2);roomDimensions(2);0];
Z = [0;0;0;0;0];
figure;
hold on;
plot3(X,Y,Z,"k","LineWidth",1.5);   % draw a square in the xy plane with z = 0
plot3(X,Y,Z+roomDimensions(3),"k","LineWidth",1.5); % draw a square in the xy plane with z = 1
set(gca,"View",[-28,35]); % set the azimuth and elevation of the plot
for k=1:length(X)-1
    plot3([X(k);X(k)],[Y(k);Y(k)],[0;roomDimensions(3)],"k","LineWidth",1.5);
end
grid on
xlabel("X (m)")
ylabel("Y (m)")
zlabel("Z (m)")
plot3(sourceCoord(1),sourceCoord(2),sourceCoord(3),"bx","LineWidth",2)
plot3(receiverCoord(1),receiverCoord(2),receiverCoord(3),"ro","LineWidth",2)
end

function h = HelperImageSource(roomDimensions,receiverCoord, ....
                    sourceCoord,A,FVect,fs,useHRTF,varargin)
% HELPERIMAGESOURCE Estimate impulse response of shoebox room
% roomDimensions: Room dimensions, specified as a row vector with three
%                 values.
% receiverCoord: Receiver coordinates, specified as a row vector with 3
%                values
% sourceCoord: Source coordinates, specified as a row vector with 3
%              values
% A:  Wall absorption coefficient matrix, specified as a L-by-6 matrix,
%     where L is the number of frequency bands.
% FVect: Vector of frequencies, of length L.
% fs: Sampling rate, in Hertz
% useHRTF: Specify as true to use HRTF interpolation
% hrtfData: Specify is useHRTF is true
%sourcePosition: Specify is useHRTF is true

hrtfData = [];
sourcePosition = [];
if useHRTF
    hrtfData = varargin{1};
    sourcePosition = varargin{2};
end

x = sourceCoord(1);
y = sourceCoord(2);
z = sourceCoord(3);
sourceXYZ = [-x -y -z; ...
             -x -y  z; ...
             -x  y -z; ...
             -x  y  z; ...
              x -y -z; ...
              x -y  z; ...
              x  y -z; ...
              x  y  z].';

Lx=roomDimensions(1); 
Ly=roomDimensions(2);
Lz=roomDimensions(3);

V = Lx*Ly*Lz;

WallXZ=Lx*Lz;
WallYZ=Ly*Lz;
WallXY=Lx*Ly;

S = WallYZ*(A(:,1)+A(:,2))+WallXZ.*(A(:,3)+A(:,4))+WallXY.*(A(:,5)+A(:,6));

c  = 343; % Speed of sound (m/s)
RT60 = (55.25/c)*V./S;

impResLength = fix(max(RT60)*fs);

impResRange=c*(1/fs)*impResLength;

nMax = min(ceil(impResRange./(2.*Lx)),10);
lMax = min(ceil(impResRange./(2.*Ly)),10);
mMax = min(ceil(impResRange./(2.*Lz)),10);

B=sqrt(1-A);

BX1=B(:,1);
BX2=B(:,2); 
BY1=B(:,3); 
BY2=B(:,4); 
BZ1=B(:,5); 
BZ2=B(:,6);

surface_coeff=[0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
q=surface_coeff(:,1).'; 
j=surface_coeff(:,2).';
k=surface_coeff(:,3).'; 

FFTLength=512; 
HalfLength=fix(FFTLength./2);
OneSidedLength = HalfLength+1;
win = hann(FFTLength+1);

FVect2=[0 FVect fs/2]';

h = zeros(impResLength,2);

for n=-nMax:nMax
    Lxn2=n*2*Lx;
    for l=-lMax:lMax
        Lyl2=l*2*Ly;

       if useHRTF
            imagesVals = zeros(FFTLength+size(hrtfData,3),2,2*lMax+1,8);
       else
           imagesVals = zeros(FFTLength+1,2,2*lMax+1,8);
       end

        Li = size(imagesVals,1);
        isDelayValid = zeros(2*lMax+1,8);
        start_index_HpV = zeros(2*lMax+1,8);
        stop_index_HpV = zeros(2*lMax+1,8);
        start_index_hV = zeros(2*lMax+1,8);

        parfor mInd=1:2*mMax+1

            m = mInd - mMax - 1;

            Lzm2=m*2*Lz;
            xyz = [Lxn2; Lyl2; Lzm2];

            isourceCoordV=xyz - sourceXYZ;
            xyzV = isourceCoordV - receiverCoord.';
            distV = sqrt(sum(xyzV.^2));
            delayV = (fs/c)*distV;

            ImagePower = BX1.^abs(n-q).*BY1.^abs(l-j).*BZ1.^abs(m-k).*BX2.^abs(n).*(BY2.^abs(l)).*(BZ2.^abs(m));
            ImagePower2 = [ImagePower(1,:); ImagePower; ImagePower(6,:)];

            ImagePower2 = ImagePower2./distV;

            validDelay = delayV<= impResLength;

            if sum(validDelay)==0
                continue;
            end

            isDelayValid(mInd,:) = validDelay;

            ImagePower2 = interp1(FVect2./(fs/2),ImagePower2,linspace(0,1,257));
            if isrow(ImagePower2)
                ImagePower2 = ImagePower2.';
            end
            ImagePower3 = [ImagePower2; conj(ImagePower2(HalfLength:-1:2,:))];

            h_ImagePower = real(ifft(ImagePower3,FFTLength));
            h_ImagePower = [h_ImagePower(OneSidedLength:FFTLength,:); h_ImagePower(1:OneSidedLength,:)];
            h_ImagePower = win.*h_ImagePower;

            if useHRTF
                hyp = sqrt(xyzV(1,:).^2+xyzV(2,:).^2);
                elevation = atan(xyzV(3,:)./(hyp+realmin));
                azimuth = atan2(xyzV(2,:),xyzV(1,:));

                desiredPosition = [azimuth.',elevation.']*180/pi;

                interpolatedIR  = interpolateHRTF(hrtfData,sourcePosition,desiredPosition,"Algorithm","VBAP");
                interpolatedIR = squeeze(permute(interpolatedIR,[3 2 1]));

                pad_ImagePower = zeros(512,2);

                for index=1:8
                    hrir0 = interpolatedIR(:,:,index);
                    hrir_ext=[hrir0; pad_ImagePower];
                    for ear=1:2
                        imagesVals(:,ear,mInd,index)=filter(h_ImagePower(:,index),1,hrir_ext(:,ear));
                    end
                end
            else
                for index=1:8
                    for ear=1:2
                        imagesVals(:,ear,mInd,index)=h_ImagePower(:,index);
                    end
                end
            end

            adjust_delay = round(delayV) - (fix(FFTLength/2))+1;

            len_h=Li;
            start_index_HpV(mInd,:) = max(adjust_delay+1+(adjust_delay>=0),1);
            stop_index_HpV(mInd,:) = min(adjust_delay+1+len_h,impResLength);
            start_index_hV(mInd,:) = max(-adjust_delay,1);

        end
        stop_index_hV = start_index_hV + (stop_index_HpV - start_index_HpV);

        for index2=1:size(imagesVals,3)
            for index3=1:8
                if isDelayValid(index2,index3)
                    h(start_index_HpV(index2,index3):stop_index_HpV(index2,index3),:)= h(start_index_HpV(index2,index3):stop_index_HpV(index2,index3),:) + squeeze(imagesVals(start_index_hV(index2,index3):stop_index_hV(index2,index3),:,index2,index3));
                end
            end
        end

    end
end

h = h./max(abs(h));
end