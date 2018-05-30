function [data,frames] = read_camera_data(path,cleanup)
%% read_camera_data - Read visual-inertial data captured by senserve
%
% Syntax:
%   [data,frames] = read_camera_data(path,cleanup)
% 
% In:
%   path    - Path to a folder with the senserve data (gz) and video file (mov)
%   cleanup - Re-compress after extraction (default: false)
% 
% Out:
%   data   - Senserve-style data matrix
%   frames - Cell-array of the extracted png frame filenames
% 
% Description:
%   Read visual(-inertial) data collected by ios senserve.
% 
% Copyright: 
%   (c) Arno Solin, 2016
% 


%% Set defaults

  % Require the path
  if nargin<1
    error('At least tha data path is required.')
  end
  
  % Remove temporary files after loading
  if nargin<2 || isempty(cleanup)
    cleanup = false;
  end

  
%% Extarct gz to csv

  % Find original gz file
  files = dir(fullfile(path,'*.gz'));

  % Return if no files found
  if length(files)<1
    error('No gz files to read in %s',path)
  end
  
  % Take first match
  filename_gz = fullfile(path,files(1).name);
  
  % Data filename
  [pathstr,name,ext] = fileparts(filename_gz);
  filename = fullfile(pathstr,name);
  
  % gunzip if file does not exist
  if ~exist(filename,'file')
    gunzip(filename_gz);
  end

  
%% Read CSV data

  % This only works if the file truly is in CSV format
  data = csvread(filename);

  % Remove the gunzipped file
  if cleanup
    delete(filename)
  end
  
%% Extract video frames

  % End here if no MOV file with matching name
  if ~exist([filename '.mov'],'file')
    return
  end

  % Only run FFMPEG if files have not been extracted
  if ~exist(fullfile(pathstr,'frame-00001.png'),'file')
  
    % Extract frames using FFMPEG (requires unix)
    status = unix(sprintf('ffmpeg -r 2 -i %s.mov %s/frame-%%05d.png', ...
      filename,pathstr));
  
    % If problems encountered, report
    if status~=0
      warning('Problem extracting frames from %s.mov',filename)
    end
  
  end

  % List files
  files = dir(fullfile(pathstr,'frame-*.png'));
  
  % Turn list into a cell array
  frames = cell(length(files),1);
  for i=1:length(files)
    frames{i} = fullfile(pathstr,files(i).name);
  end
  
  % Remove the extracted frames
  if cleanup
    delete(fullfile(pathstr,'frame-*.png'))
  end
  
  
%% Sanity check that we have the timings

  % The frame indices
  ind = (data(:,2)==22);
  
  % Throw error if the frame count does not match
  if sum(ind)~=(length(frames)-1)
    error('Video frame count %i-1 does not match the timing count %i.', ...
        length(frames),sum(ind))
  end
  

  
  
  