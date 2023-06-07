            % Fuzzy C-means classification with 2 classes
              dsFrame = readrawfile(fullfile(FolderDir, imageLists{f}));
              imwrite(dsFrame,[DataDir, 'Evaluation\temp.tiff'])

            data = (movingObject(:)); % data array
            [center,U,obj_fcn] = fcm(data,7); 

            Finding the pixels for each class
            maxU = max(U);
            index1 = find(U(1,:) == maxU);
            index2 = find(U(2,:) == maxU);
            index3 = find(U(3,:) == maxU);
            index4 = find(U(4,:) == maxU);
            index5 = find(U(5,:) == maxU);
            index6 = find(U(6,:) == maxU);
            index7 = find(U(7,:) == maxU);

            Assigning pixel to each class by giving them a specific value
            fcmImage(1:length(data))=0.3;       
            fcmImage(index1)= 0.0;
            fcmImage(index2)= 1;
            fcmImage(index3)= 0.0;
            fcmImage(index4)= 0.4;
            fcmImage(index5)= 0.2;
            fcmImage(index6)= 0.1;
            fcmImage(index7)= 0.0;

            Reshapeing the array to a image
            imagNew = reshape(fcmImage,1440,1560);
            figure(101);imshow(imagNew,[]);impixelinfo; 