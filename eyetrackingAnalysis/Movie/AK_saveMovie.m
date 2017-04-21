function AK_saveMovie( filename, Mstruct, frameRate, quality )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% check inputs
if nargin < 2
    error('AK_saveMovie requires at least two inputs: the fullfile directory (string); and an array of movie frames (struct, use getframes)')
end
if nargin < 3
    frameRate = [];
end
if nargin < 4
    quality = [];
end

% save movie w/ settings
myMovie = VideoWriter(filename);
if ~isempty(frameRate)
    myMovie.FrameRate = frameRate; % Default 30
end
if ~isempty(quality)
    myMovie.Quality = quality; % Default 75
end
% write video
open(myMovie);
writeVideo(myMovie,Mstruct);
close(myMovie);

% movie2avi(M,fullfile(data_dir,[edfName '.avi']),'fps',60)

end

