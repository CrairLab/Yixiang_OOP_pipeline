function varargout = GUI_dimReduction(varargin)
% GUI_DIMREDUCTION MATLAB code for GUI_dimReduction.fig
%      GUI_DIMREDUCTION, by itself, creates a new GUI_DIMREDUCTION or raises the existing
%      singleton*.
%
%      H = GUI_DIMREDUCTION returns the handle to a new GUI_DIMREDUCTION or the handle to
%      the existing singleton*.
%
%      GUI_DIMREDUCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DIMREDUCTION.M with the given input arguments.
%
%      GUI_DIMREDUCTION('Property','Value',...) creates a new GUI_DIMREDUCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_dimReduction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_dimReduction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Author: Yixiang Wang
%      Contact: yixiang.wang@yale.edu
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_dimReduction

% Last Modified by GUIDE v2.5 18-Jun-2019 20:21:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_dimReduction_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_dimReduction_OutputFcn, ...
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


% --- Executes just before GUI_dimReduction is made visible.
function GUI_dimReduction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_dimReduction (see VARARGIN)

% Choose default command line output for GUI_dimReduction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_dimReduction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_dimReduction_OutputFcn(hObject, eventdata, handles) 
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
% Load movie button.

    %load dFoF movie using uiopen
    set(handles.edit5, 'Visible', 'On')
    set(handles.edit5, 'String', 'Loading...')
    
    uiopen('Please load the 3D dF over F movie!')
    
    %Store the movie into a variable curMovie
    try
        vList = whos; 
        for i = 1:size(vList,1)
            %Search for 3D matrix
            if length(vList(i).size) == 3 
                curMovie = eval(vList(i).name);
                hObject.UserData = curMovie;
                set(handles.edit5, 'String', 'Finished!')
                break
            end
        end  
    catch
        set(handles.edit5, 'String', 'Error!')
        warning('Can not load the dF over F movie!')        
    end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Popupmenu to choose whether do pixelwise or framewise analysis

%Load existed initial parameters
iniParameters = get(handles.pushbutton2, 'UserData');
contents = cellstr(get(hObject, 'String'));
curString = contents{get(hObject,'Value')};
if strcmp(curString, 'Pixelwise')
        iniParameters.tflag = 0;
elseif strcmp(curString, 'Framewise')
        iniParameters.tflag = 1;
end

%Update initial parameters
set(handles.pushbutton2, 'UserData', iniParameters);

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


% --- Executes on key press with focus on popupmenu1 and none of its controls.
function popupmenu1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pushbutton to load an existed dimReduction object

set(handles.edit13, 'Visible', 'On')
set(handles.edit13, 'String', 'Wait')

%Open a saved dimReduction object using dialog window
try
    uiopen('Please load a dimReduction object');
catch
    warning('Please load a dimReduction object')
end

%Store the loaded object to a variable curObj
vList = whos; 
for i = 1:size(vList,1)
    if strcmp(vList(i).class, 'dimReduction')
        curObj = eval(vList(i).name);
        set(handles.output, 'UserData', curObj)
        set(handles.edit13, 'Visible', 'Off')
        break
    end
end
set(handles.output, 'UserData', curObj);

%Redisplay parameters;
displayParam(curObj, handles);

%Plot tSNE and diffusion map
plotAxes(curObj, handles);

% --- Executes on selection change in popupmenu1.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pushbutton to load default parameters

%Set default parameters
%pixelwise analysis
iniParameters.tflag = 0;
set(handles.popupmenu1,'Value',2);
%Do not constraint on location
iniParameters.locflag = 0;
set(handles.checkbox1,'Value',0);
iniParameters.locfactor = 0;
set(handles.edit1,'Value',0);
%Downsample by 2 both temporally and spatially
iniParameters.fd = [2, 2];
set(handles.edit2,'Value',2);
set(handles.edit2,'String',num2str(2));
set(handles.edit3,'Value',2);
set(handles.edit3,'String',num2str(2));
%Save default parameters;
hObject.UserData = iniParameters;



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Checkbox to choose whether inject location information into analysis

% Hint: get(hObject,'Value') returns toggle state of checkbox1
iniParameters = get(handles.pushbutton2, 'UserData');
curValue = get(hObject, 'Value');
if curValue == 0 
        iniParameters.locflag = 0;
elseif curValue == 1
        iniParameters.locflag = 1;
end
set(handles.pushbutton2, 'UserData', iniParameters);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
iniParameters = get(handles.pushbutton2, 'UserData');
iniParameters.locfactor = str2double(get(hObject, 'String'));
set(handles.pushbutton2, 'UserData', iniParameters);

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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pushbutton to run dimReduction analysis/construct dimReduction object

%Show progress
set(handles.edit4, 'Visible', 'On')
set(handles.edit4, 'String', 'Running...')

%Load parameters
param = get(handles.pushbutton2, 'UserData');
curMovie = get(handles.pushbutton1, 'UserData');

%Run dimReduction
try
    curObj = dimReduction(curMovie, param.tflag, param.locflag,...
        param.locfactor, param.fd);
    set(handles.edit4, 'String', 'Finished!')
catch
    set(handles.edit4, 'String', 'Error!')
end
set(handles.output, 'UserData', curObj);

%Redisplay parameters;
displayParam(curObj, handles)

%Plot tSNE and diffusion map
plotAxes(curObj, handles);

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
iniParameters = get(handles.pushbutton2, 'UserData');
iniParameters.fd(1) = str2double(get(hObject, 'String'));
set(handles.pushbutton2, 'UserData', iniParameters);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
iniParameters = get(handles.pushbutton2, 'UserData');
iniParameters.fd(2) = str2double(get(hObject, 'String'));
set(handles.pushbutton2, 'UserData', iniParameters);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Pushbutton to load default parameters
iniParameters.tflag = 0;
%Do not constraint on location
iniParameters.locflag = 0;
iniParameters.locfactor = 0;
%Downsample by 2 both temporally and spatially
iniParameters.fd = [2, 2];
%Save default parameters;
hObject.UserData = iniParameters;



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pushbutton to renew tSNE analysis using new parameters

%Show progress
set(handles.edit8, 'Visible', 'On')
set(handles.edit8, 'String', 'Running...')

%Load renewed parameters
tParam = get(handles.uipanel2, 'UserData');
tParam.px = str2double(get(handles.edit6, 'String'));
tParam.Exaggeration = str2double(get(handles.edit7, 'String'));
contents = cellstr(get(handles.popupmenu2, 'String'));
curString = contents{get(handles.popupmenu2,'Value')};
tParam.Distance = curString;
set(hObject.Parent, 'UserData', tParam);

%Redo tSNE using new parameters
curObj = get(handles.output, 'UserData');
curObj.tParam = tParam;
curObj.Y = dimReduction.doTSNE(curObj.A_rd, tParam);
set(handles.output, 'UserData', curObj);
set(handles.edit8, 'String', 'Finished!')

%Redisplay parameters;
displayParam(curObj, handles);

%Renew the current plots
plotAxes(curObj, handles);


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Load default tSNE parameters
tParam.Distance = 'euclidean';
tParam.nd = 3;
tParam.np = 50;
tParam.px = 30;
tParam.stdFlag = false;
tParam.vflag = 1;
tParam.Exaggeration = 8;
tParam.options = statset('MaxIter', 500);
hObject.UserData = tParam;



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pushbutton to renew diffusion map analysis using new parameters

%Show progress
set(handles.edit12, 'Visible', 'On')
set(handles.edit12, 'String', 'Running...')

%Get renewed parameters
dParam = get(handles.uipanel3, 'UserData');
dParam.t = str2double(get(handles.edit9, 'String'));
dParam.m = str2double(get(handles.edit10, 'String'));
dParam.sigma = str2double(get(handles.edit11, 'String'));
set(hObject.Parent, 'UserData', dParam);

%Redo diffusion map using new parameters
curObj = get(handles.output, 'UserData');
[curObj.Dmap, curObj.dParam, ~] = ...
    dimReduction.diffmap(curObj.A_rd, dParam.t, dParam.m, dParam.sigma);
set(handles.output, 'UserData', curObj);
set(handles.edit12, 'String', 'Finished!')

%Redisplay parameters;
displayParam(curObj, handles);

%Renew the current plots
plotAxes(curObj, handles);


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Load default parameters for diffusion maps
dParam.t = 2;
dParam.m = 3;
dParam.sigma = [];
hObject.UserData = dParam;



function displayParam(curObj, handles)
% display parameters in edit-text boxes
% Inputs:
% curObj     input dimReduction object
% handles    handles of UI components

%Display parameters for dimReduction
set(handles.popupmenu1,'Value', curObj.tflag + 2);
%Do not constraint on location
set(handles.checkbox1,'Value',curObj.locflag);
set(handles.edit1,'Value',0);
%Downsample by 2 both temporally and spatially
set(handles.edit2,'String',num2str(curObj.fd(1)));
set(handles.edit3,'String',num2str(curObj.fd(2)));

%Display tSNE parameters
tParam = curObj.tParam;
set(handles.edit6, 'String', num2str(tParam.px));
set(handles.edit7, 'String', num2str(tParam.Exaggeration));

contents = cellstr(get(handles.popupmenu2,'String'));
curVal = find(strcmp(contents, tParam.Distance));
set(handles.popupmenu2, 'Value', curVal);

%Display diffusion map parameters
dParam = curObj.dParam;
set(handles.edit9, 'String', num2str(dParam.t));
set(handles.edit10, 'String', num2str(dParam.m));
set(handles.edit11, 'String', num2str(dParam.sigma));




function plotAxes(curObj, handles)
% plot in Axes (renew tSNE map and diffusion map)
% Inputs:
% curObj     input dimReduction object
% handles    handles of UI components

curY = curObj.Y;
curDmap = curObj.Dmap;

cmap = 1:size(curY,1);
scatter3(handles.axes1, curY(:,1),curY(:,2),curY(:,3),[],cmap,'filled');

cmap = 1:size(curDmap,1);
scatter3(handles.axes2, curDmap(:,1),curDmap(:,2),curDmap(:,3),[],cmap,'filled');


function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curObj = get(handles.output, 'UserData');
CorrespondMaps(curObj);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pushbutton to save current dimReduction to a .mat files
c = clock;
timetag = ['_' num2str(c(1)) num2str(c(2)) num2str(c(3)) num2str(c(4)) num2str(c(5))];
curObj = get(handles.output, 'UserData');
paramtag = ['_' num2str(curObj.tflag) '_' num2str(curObj.locflag), '_', ...
    num2str(curObj.locfactor), '_', num2str(curObj.fd(1)) num2str(curObj.fd(2))];

if curObj.tflag
    methodtag = '_Framewise_';
else
    methodtag = '_Pixelwise_';
end

uisave('curObj', ['dimReduction' methodtag paramtag timetag '.mat'])
