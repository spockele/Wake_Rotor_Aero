clear; 
%% UNSTEADY
folderName = './figures/field/unsteady/K_0.10';
% load the images
folder = dir([folderName '/*.jpg']);
N = size(folder,1);
images = cell(N,1);
for ii=1:N; images{ii} = imread(strcat(folderName,sprintf('/%.0f.jpg',ii)));    end

% %% STEADY
% folderName = './figures/field/steady';
% load(strcat(folderName,'/subplot_figList.mat'));
% N = size(textList,1);
% images = cell(N,1);
% for ii=1:N; images{ii} = imread(textList{ii});    end

%% create the video writer with # fps
 writerObj = VideoWriter(strcat(folderName,'video.avi'));
 writerObj.FrameRate = 2;
 % set the seconds per image
 secsPerImage = ones(1,N);
 % open the video writer
 open(writerObj);
 % write the frames to the video
 for u=1:N
     % convert the image to a frame
     frame = im2frame(images{u});
     for v=1:secsPerImage(u) 
         writeVideo(writerObj, frame);
     end
 end
 % close the writer object
 close(writerObj);