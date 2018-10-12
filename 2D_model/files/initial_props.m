function varargout = initial_props(varargin)
% INITIAL_PROPS M-file for initial_props.fig
%      INITIAL_PROPS, by itself, creates a new INITIAL_PROPS or raises the existing
%      singleton*.
%
%      H = INITIAL_PROPS returns the handle to a new INITIAL_PROPS or the handle to
%      the existing singleton*.
%
%      INITIAL_PROPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INITIAL_PROPS.M with the given input arguments.
%
%      INITIAL_PROPS('Property','Value',...) creates a new INITIAL_PROPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before initial_props_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to initial_props_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help initial_props

% Last Modified by GUIDE v2.5 13-Jun-2009 18:38:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @initial_props_OpeningFcn, ...
                   'gui_OutputFcn',  @initial_props_OutputFcn, ...
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


% --- Executes just before initial_props is made visible.
function initial_props_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to initial_props (see VARARGIN)

% Choose default command line output for initial_props
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes initial_props wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = initial_props_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function lx_Callback(hObject, eventdata, handles)
% hObject    handle to lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lx as text
%        str2double(get(hObject,'String')) returns contents of lx as a double


% --- Executes during object creation, after setting all properties.
function lx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ly_Callback(hObject, eventdata, handles)
% hObject    handle to ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ly as text
%        str2double(get(hObject,'String')) returns contents of ly as a double


% --- Executes during object creation, after setting all properties.
function ly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_Callback(hObject, eventdata, handles)
% hObject    handle to n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n as text
%        str2double(get(hObject,'String')) returns contents of n as a double


% --- Executes during object creation, after setting all properties.
function n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rad_Callback(hObject, eventdata, handles)
% hObject    handle to rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rad as text
%        str2double(get(hObject,'String')) returns contents of rad as a double


% --- Executes during object creation, after setting all properties.
function rad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m_Callback(hObject, eventdata, handles)
% hObject    handle to m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m as text
%        str2double(get(hObject,'String')) returns contents of m as a double


% --- Executes during object creation, after setting all properties.
function m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rm.
function rm_Callback(hObject, eventdata, handles)
% hObject    handle to rm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rm


% --- Executes on button press in rr.
function rr_Callback(hObject, eventdata, handles)
if get(handles.rr,'value')
    set(handles.rpm,'value',false);
else
    set(handles.rpm,'value',true);
end
    

% --- Executes on button press in rpm.
function rpm_Callback(hObject, eventdata, handles)
if get(handles.rpm,'value')
    set(handles.rr,'value',false);
else
    set(handles.rr,'value',true);
end


% --- Executes on button press in rc.
function rc_Callback(hObject, eventdata, handles)
% hObject    handle to rc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rc


% --- Executes on button press in hr.
function hr_Callback(hObject, eventdata, handles)
% hObject    handle to hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hr


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);


