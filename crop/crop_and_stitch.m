function varargout = crop_and_stitch(varargin)
% CROP_AND_STITCH MATLAB code for crop_and_stitch.fig
%      CROP_AND_STITCH, by itself, creates a new CROP_AND_STITCH or raises the existing
%      singleton*.
%
%      H = CROP_AND_STITCH returns the handle to a new CROP_AND_STITCH or the handle to
%      the existing singleton*.
%
%      CROP_AND_STITCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CROP_AND_STITCH.M with the given input arguments.
%
%      CROP_AND_STITCH('Property','Value',...) creates a new CROP_AND_STITCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before crop_and_stitch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to crop_and_stitch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help crop_and_stitch

% Last Modified by GUIDE v2.5 20-Jun-2019 08:53:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @crop_and_stitch_OpeningFcn, ...
                   'gui_OutputFcn',  @crop_and_stitch_OutputFcn, ...8
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before crop_and_stitch is made visible.
function crop_and_stitch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crop_and_stitch (see VARARGIN)


if exist('preferences.mat', 'file')
    load preferences.mat
end
if exist('crop_size', 'var')
    set(handles.ed_blocky, 'String', num2str(crop_size(1)));
    set(handles.ed_blockx, 'String', num2str(crop_size(2)));
    set(handles.ed_blockz, 'String', num2str(crop_size(3)));
end

if exist('grid_dim', 'var')
    set(handles.ed_gridy, 'String', num2str(grid_dim(1)));
    set(handles.ed_gridx, 'String', num2str(grid_dim(2)));
    set(handles.ed_gridz, 'String', num2str(grid_dim(3)));
end

if exist('overlap', 'var')
    set(handles.ed_overlapy, 'String', num2str(overlap(1)));
    set(handles.ed_overlapx, 'String', num2str(overlap(2)));
    set(handles.ed_overlapz, 'String', num2str(overlap(3)));
    
    set(handles.ed_soverlapy, 'String', num2str(overlap(1)));
    set(handles.ed_soverlapx, 'String', num2str(overlap(2)));
    set(handles.ed_soverlapz, 'String', num2str(overlap(3)));
end

% Choose default command line output for crop_and_stitch
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crop_and_stitch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = crop_and_stitch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ed_blockx_Callback(hObject, eventdata, handles)
% hObject    handle to ed_blockx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_blockx as text
%        str2double(get(hObject,'String')) returns contents of ed_blockx as a double

function ed_blocky_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ed_blockx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_blockx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_blockz_Callback(hObject, eventdata, handles)

function ed_overlapy_Callback(hObject, eventdata, handles)
overlap = get(hObject, 'String');
set(handles.ed_overlapx, 'String', overlap);
set(handles.ed_overlapz, 'String', overlap);

function ed_overlapx_Callback(hObject, eventdata, handles)
overlap = get(hObject, 'String');
set(handles.ed_overlapz, 'String', overlap);

function ed_overlapz_Callback(hObject, eventdata, handles)

% --- Executes on button press in rb_saveall.
function rb_saveall_Callback(hObject, eventdata, handles)
% hObject    handle to rb_saveall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_saveall
if get(hObject,'Value') == 0
    set(handles.ed_sum, 'Enable', 'on')
    set(handles.ed_var, 'Enable', 'on')
else
    set(handles.ed_sum, 'Enable', 'off')
    set(handles.ed_var, 'Enable', 'off')
end

function rb_edge_discard_Callback(hObject, eventdata, handles)
function rb_edge_padding_Callback(hObject, eventdata, handles)
function rb_edge_cover_all_Callback(hObject, eventdata, handles)

function ed_sum_Callback(hObject, eventdata, handles)

function ed_var_Callback(hObject, eventdata, handles)

function ed_gridx_Callback(hObject, eventdata, handles)

function ed_gridy_Callback(hObject, eventdata, handles)

function ed_gridz_Callback(hObject, eventdata, handles)

function ed_soverlapy_Callback(hObject, eventdata, handles)
overlap = get(hObject, 'String');
set(handles.ed_soverlapx, 'String', overlap);
set(handles.ed_soverlapz, 'String', overlap);

function ed_soverlapx_Callback(hObject, eventdata, handles)
overlap = get(hObject, 'String');
set(handles.ed_soverlapz, 'String', overlap);
function ed_soverlapz_Callback(hObject, eventdata, handles)


% --- Executes on button press in cb_delete_src.
function cb_delete_src_Callback(hObject, eventdata, handles)
% hObject    handle to cb_delete_src (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_delete_src


% --- Executes on button press in btn_stitch.
function btn_stitch_Callback(hObject, eventdata, handles)
% hObject    handle to btn_stitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gx = floor(str2double(get(handles.ed_gridx, 'String')));
gy = floor(str2double(get(handles.ed_gridy, 'String')));
gz = floor(str2double(get(handles.ed_gridz, 'String')));

olp_x = str2double(get(handles.ed_soverlapx, 'String'));
olp_y = str2double(get(handles.ed_soverlapy, 'String'));
olp_z = str2double(get(handles.ed_soverlapz, 'String'));

grid = [gy, gx, gz];
overlap = [olp_y, olp_x, olp_z];
delete_source = get(handles.cb_delete_src,'Value');
stitching3d_fn(grid, overlap, delete_source);





% --- Executes on button press in btn_crop.
function btn_crop_Callback(hObject, eventdata, handles)
% hObject    handle to btn_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
by = floor(str2double(get(handles.ed_blocky, 'String')));
bx = floor(str2double(get(handles.ed_blockx, 'String')));
bz = floor(str2double(get(handles.ed_blockz, 'String')));

oy = str2double(get(handles.ed_overlapy, 'String'));
ox = str2double(get(handles.ed_overlapx, 'String'));
oz = str2double(get(handles.ed_overlapz, 'String'));

edge_mode = 0;
if get(handles.rb_edge_discard, 'Value')
    edge_mode = 1;
elseif get(handles.rb_edge_padding, 'Value')    edge_mode = 2;
elseif get(handles.rb_edge_cover_all, 'Value')
    edge_mode = 3;
end
    
save_all = get(handles.rb_saveall, 'Value');
if ~save_all
    pix_sum = str2double(get(handles.ed_sum, 'String'));
    pix_var = str2double(get(handles.ed_var, 'String'));
    crop3d_fn([by, bx, bz], [oy, ox, oz], save_all, edge_mode, pix_sum, pix_var);
    
else
    grid_dim = crop3d_fn([by, bx, bz], [oy, ox, oz], save_all, edge_mode);
end
