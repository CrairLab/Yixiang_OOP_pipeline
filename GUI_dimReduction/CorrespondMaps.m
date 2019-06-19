function varargout = CorrespondMaps(varargin)
% CORRESPONDMAPS MATLAB code for CorrespondMaps.fig
%      CORRESPONDMAPS, by itself, creates a new CORRESPONDMAPS or raises the existing
%      singleton*.
%
%      H = CORRESPONDMAPS returns the handle to a new CORRESPONDMAPS or the handle to
%      the existing singleton*.
%
%      CORRESPONDMAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRESPONDMAPS.M with the given input arguments.
%
%      CORRESPONDMAPS('Property','Value',...) creates a new CORRESPONDMAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorrespondMaps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorrespondMaps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Author: Yixiang Wang
%      Contact: yixiang.wang@yale.edu
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorrespondMaps

% Last Modified by GUIDE v2.5 18-Jun-2019 21:50:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrespondMaps_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrespondMaps_OutputFcn, ...
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


% --- Executes just before CorrespondMaps is made visible.
function CorrespondMaps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorrespondMaps (see VARARGIN)

% Choose default command line output for CorrespondMaps
if size(varargin,1) > 0
    try
        curObj = varargin{1};
        hObject.UserData = curObj;
    catch
        warning('Please load dimReuction object!')
    end
end

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CorrespondMaps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CorrespondMaps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(get(handles.output, 'UserData'))
    try
        uiopen('Please load a dimReduction object');
        set(handles.edit1, 'String', 'Object loaded from dialog window!')        
    catch
        warning('Please load a dimReduction object')
    end
    
    vList = whos; 
    for i = 1:size(vList,1)
        if strcmp(vList(i).class, 'dimReduction')
            curObj = eval(vList(i).name);
            set(handles.output, 'UserData', curObj);
            break
        end
    end
else
    set(handles.edit1, 'String', 'Object loaded from previous GUI!')   
end
set(handles.edit1, 'ForegroundColor', 'Red' )

function runCorrespondMaps(curObj, handles, methodflag, mapflag)
% Correspond dimReduction map with reference brain map or frame(s)
% Inputs:
%   curObj    input dimReduction object
%   handles   UI component handles
%   methodflag     0 for pixelwise, 1 for framewise
%   mapflag   1 for tSNE, 2 for diffusion map
if nargin<3
    mapflag = 1;
end

A_mean = curObj.A_ref;
imshow(mat2gray(A_mean), 'Parent', handles.axes1);
hold(handles.axes1, 'on');

if mapflag == 1
    Embedding = curObj.Y;
elseif mapflag == 2
    Embedding = curObj.Dmap;
end


cmap = 1:size(Embedding,1);
fig1 = figure;
map_handle = scatter3(Embedding(:,1),Embedding(:,2),Embedding(:,3),[],cmap,'filled');
colorbar;
%map_handle = plot3(Embedding(:,1),Embedding(:,2),Embedding(:,3));

data.map_handle = map_handle;
data.handles = handles;
fig1.UserData =  data;

if ~methodflag
    h = brush(fig1);
    set(h, 'Enable', 'On', 'ActionPostCallback', @callbackGetSelectedData);
elseif methodflag
    h = brush(fig1);
    set(h, 'Enable', 'On', 'ActionPostCallback', @callbackClickA3DPoint);
end
    
function callbackClickA3DPoint(src, ~)
    disp('Please only select one point!')
    data = src.UserData;
    brush_data = get(data.map_handle, 'BrushData');
    brushed_idx = find(brush_data);
    firstIdx = brushed_idx(1);
    data.brushed = firstIdx;
    src.UserData = data;
    handles = data.handles;
    try
        curObj = get(handles.output, 'UserData');
        sz_fd = curObj.sz_fd;
        A_rcs = nan(sz_fd(1)*sz_fd(2), sz_fd(3));
        A_rd = curObj.A_rd;
        
        if curObj.tflag == 1
            A_rd = A_rd';
        end
        
        try
            A_rcs(curObj.subIdx,:) = A_rd;
        catch
            A_rcs(curObj.subIdx,:) = A_rd';
        end
        
        A_rcs = reshape(A_rcs, sz_fd);
        correspondFrame = A_rcs(:,:,firstIdx);
        imshow(mat2gray(correspondFrame), 'Parent', handles.axes1);
        set(handles.edit1, 'String', ['Frame #' num2str(firstIdx)])
    catch
        warning('Can not display corresponding frame!')
    end
    


function callbackGetSelectedData(src, ~)
    data = src.UserData;
    brush_data = get(data.map_handle, 'BrushData');
    brushed_idx = find(brush_data);
    data.brushed = brushed_idx;
    src.UserData = data;
    handles = data.handles;
    try
        curObj = get(handles.output, 'UserData');
        xy_sub = curObj.xy_sub;
        %imshow(mat2gray(A_mean), 'Parent', handles.axes1);
        A_mean = curObj.A_ref;
        imshow(mat2gray(A_mean), 'Parent', handles.axes1);
        hold(handles.axes1, 'on');
        plot(handles.axes1,xy_sub(brushed_idx,2),xy_sub(brushed_idx,1),'r.')
        hold(handles.axes1, 'off');
    catch
        warning('Can not load xy subscripts!')
    end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton1.
function pushbutton1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)






function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curObj = get(handles.output, 'UserData');
tflag = curObj.tflag;
mapflag = get(handles.popupmenu1, 'Value');

if tflag == 0 %Pixelwise
    set(handles.edit1, 'String', 'Pixelwise correspondence!')
elseif tflag == 1 %Framewise
    set(handles.edit1, 'String', 'Framewise correspondence!')
end

runCorrespondMaps(curObj, handles, tflag, mapflag)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
