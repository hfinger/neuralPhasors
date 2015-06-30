
function varargout = SimulationGUI(varargin)
% SIMULATIONGUI MATLAB code for SimulationGUI.fig
%      SIMULATIONGUI, by itself, creates a new SIMULATIONGUI or raises the existing
%      singleton*.
%
%      H = SIMULATIONGUI returns the handle to a new SIMULATIONGUI or the handle to
%      the existing singleton*.
%
%      SIMULATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMULATIONGUI.M with the given input arguments.
%
%      SIMULATIONGUI('Property','Value',...) creates a new SIMULATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SimulationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SimulationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SimulationGUI

% Last Modified by GUIDE v2.5 04-Jun-2015 18:00:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SimulationGUI_OpeningFcn, ...
    'gui_OutputFcn',  @SimulationGUI_OutputFcn, ...
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

% --- configuration function
function param = config()
param.class_path = '/net/store/nbp/projects/phasesim/src_kstandvoss/matlab/classes';
param.work_path = '/net/store/nbp/projects/phasesim/workdir';

% --- Executes just before SimulationGUI is made visible.
function SimulationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SimulationGUI (see VARARGIN)

% number of added jobs
global num_jobs;
num_jobs = 0;
% currently selected job
global selected;
selected = 0;
% status refreshing
global running;
running = false;

p = config;
% Choose default command line output for SimulationGUI
handles.output = hObject;

% populate popupmenu with choosable gridjobs
allFiles = dir(p.class_path);
allNames = {allFiles(~[allFiles.isdir]).name};
set(handles.popupmenu1,'String',allNames);

% populate table with respective parameters and values
list = get(handles.popupmenu1,'String');
val = list{get(handles.popupmenu1,'Value')};
val = val(1:end-2); %selected Gridjobclass
make = str2func(val);
obj = make(); %create object
[fields, values] = build_table(hObject,handles,obj,val);

set(handles.table1,'data',[fields, values']);
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes SimulationGUI wait for user response (see UIRESUME)
%uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = SimulationGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.table1.Position(3) = handles.table1.Extent(3);
handles.status.Position(3) = handles.status.Extent(3);
set(handles.Add_job,'tooltipString','Add job to modify parameters');
set(handles.load,'tooltipString','Load job from Disc.mat');
set(handles.status,'tooltipString','Leftclick on Error/Log file to open');

varargout{1} = handles.output;

% --- Executes on selection change in popupmenu1 to select other job.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

list = get(handles.popupmenu1,'String');
val = list{get(handles.popupmenu1,'Value')};
val = val(1:end-2); %new selected job
make = str2func(val);
obj = make();
[fields values] = build_table(hObject,handles,obj,val);
set(handles.table1,'data',[fields, values']);
set(handles.table1,'ColumnEditable',[false false false]);


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



% --- Executes on button press in 'run' button.
% --- Starts execution of added gridjobs
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global allJobs; %cell-array of added jobs
global running;

for k=1:size(allJobs,2)
    allParams{k} = allJobs{k}.params;
end

if size(allJobs,2) > 0;
    clc;
    gridjob = Gridjob(allParams);
    start(gridjob); %use start-function of Gridjobclass to begin execution
    if ~running && ~gridjob.params.Gridjob.runLocal
        show_status(hObject,handles,allJobs); %show status of execution
    end
end

% --- Updates Status table every 5 seconds to show job status
function show_status(hObject,handles,jobs)

[~,user] = system('whoami'); %find executing user
con = config;

global constructedFromFolder; %workdir
global running; %variable to disable actualization

newsize = 0;
oldsize = 0; %only ckeck for finished jobs if status size changed 

running = true;
while running
    [~,status] = system('qstat'); %get status
    
    newsize = size(status);
    myids = cell(1);
    mystates = cell(1);
    mynames = cell(1);
    tasks = cell(1);
    queues = cell(1);
    errorfiles = cell(1);
    logfiles = cell(1);
    
    if ~isempty(status)
        C = textscan(status,'%s %d %s %s %s %s %s %s %s %s','HeaderLines',2);
        jobids = C{1};
        states = C{5};
        names = C{3};
        %-- find only jobs submitted by user
        k = 1;
        for i=1:size(states,1)
            if ~isempty(strfind(user,C{4}{i}))
                myids{k} = jobids{i};
                mystates{k} = states{i};
                mynames{k} = names{i};
                %-- find error and logfiles
                if strcmp(states{i},'r') || strcmp(states{i},'e')
                    tasks{i} = C{10}{k};
                    queues{i} = C{8}{k};
                    %name of errofile
                    filename = strcat('temp_',char(C{3}{i}),'/',char(C{3}{i}),'.e',char(myids{k}),'.',tasks{i});
                    %name of logfile
                    logname = strcat('temp_',char(C{3}{i}),'/',char(C{3}{i}),'.o',char(myids{k}),'.',tasks{i});
                    dat = java.io.File(filename);
                    if length(dat) >0
                        mystates{k} = 'e';
                    end
                    errorfiles{i} = filename;
                    logfiles{i} = logname;
                else
                    tasks{i} = C{9}{k};
                    queues{i} = '-';
                    errorfiles{i} = '-';
                    logfiles{i} = '-';
                end
                k = k+1;
            end
        end
    else
        running = false;
    end
    
    if newsize < oldsize
        %--find finished jobs
        for j=1:size(jobs,2)
            job = jobs{j};
            name = job.params.Gridjob.jobname;
            filename = fullfile(con.work_path,constructedFromFolder,strcat('temp_',char(name)));
            allFiles = dir(filename{1});
            allNames = {allFiles.name};
            %look for file isfinished - if its deleted the job has finished
            if ~ismember(allNames,'isfinished')
                count_string = strfind(allNames,'.o'); %get number of tasks
                occur = 0;
                nums = cell(1);
                for l = 1:size(count_string,2)
                    if ~isempty(count_string{l})
                        occur = occur+1;
                        nums = regexp(allNames{l},'o\d+','match'); %get jobid
                    end
                end
                if ~isempty(nums)
                    nums = nums{1};
                    nums = nums(2:end);
                    myids = [nums,myids];
                    mynames = [name,mynames];
                    mystates = ['f',mystates];
                    tasks = [strcat('1:',num2str(occur)),tasks];
                    queues = ['-',queues];
                    errorfiles = ['-',errorfiles];
                    logfiles = ['-',logfiles];
                end
            end
        end
    end
    set(handles.status,'data',[myids',mynames',mystates',tasks',queues',errorfiles',logfiles']);
    oldsize = newsize;
    pause(5);
end


% --- Populates Table with chosen parameters and respective values
% ob       class object to fill table with
% val      name of the chosen subclass
function [fields, values] = build_table(hobject,handles,obj,val)

gridfields = fieldnames(obj.params.Gridjob); % gridjob-superclass parameters
subfields = fieldnames(obj.params.(val));% subclass parameters
[subfields,subvalues] = set_values(subfields,val,false,obj);
[gridfields,gridvalues] = set_values(gridfields,'Gridjob',false,obj);

fields = [subfields;gridfields]; % parameterlist
values = [subvalues gridvalues];

% --- Gets parameter-values for populating table
function [fields,values] = set_values(fields,str,struct,obj)

%find values for each field
k = 1;
while k <= size(fields,1)
    string = fields(k); %get fieldname
    %--if superfield was struct get structfield parameters
    if struct
        result = strcat('obj.params.',str,'.',string{1});
        result = strcat('elem = ',result,';');
        eval(result);
        %--else just take fieldvalue
    else
        elem = obj.params.(str).(string{1});
    end
    %-- if no standard value is set
    if isempty(elem)
        elem = '-';
    end
    %-- if value is again a struct call function recursivly to fill struct
    %fields
    if isstruct(elem)
        structfields = fieldnames(elem)
        if ~isempty(structfields)
            %-- recursive call for struct values
            [structfields,structvalues] = set_values(structfields,obj,strcat(str,'.',string{1}),true,obj);
            structfields = strcat(str,'.',string{1},'.',structfields); %add struct name to fieldnames
            fields = vertcat(fields{1:k},structfields,fields{k+1:end}); %add structfields to other fieldnames
            fields(k) = strcat(str,'.',fields(k));
            
            elem = 'struct'; %struct is just denoted as such
            values{k} = elem;
            values= horzcat(values{1:k},structvalues,values{k+1:end}); %add structvalues to values
            %-- jump to next field after structfields
            k = k+size(structfields,1)+1;
            continue;
        else
            elem = 'struct';
        end
        % -- write logical as words instead of 0 or 1
    elseif islogical(elem)
        elem_string = '';
        for j=1:size(elem,2)
            if elem(j)
                elem_string = strcat(elem_string,'true ');
            else
                elem_string = strcat(elem_string,'false ');
            end
        end
        elem = elem_string;
        % -- write numbers as strings
    elseif isnumeric(elem)
        if size(elem,2)>0
            elem = num2str(elem);
        end
        % -- write cell arrays as strings (else can't be displayed in uitable)
    elseif iscell(elem)
        result = '{';
        for j=1:size(elem,2)
            result = strcat(result,num2str(elem{j}),',');
        end
        elem = strcat(result(1:end-1),'}');
    end
    
    values{k} = elem;   % get parametervalues of gridjobclass
    
    if ~struct
        fields(k) = strcat(str,'.',fields(k));
    end
    k = k +1;
end

% --- Executes when entered data in editable cell(s) in table1.
% --- Modifiable table with job-parameters
function table1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to table1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

global allJobs;
global selected;

list = get(handles.popupmenu1,'String');
val = list{get(handles.popupmenu1,'Value')};
val = val(1:end-2);
gridfields = fieldnames(allJobs{selected}.params.Gridjob); % gridjob-superclass parameters
subfields = fieldnames(allJobs{selected}.params.(val)); % subclass parameters
fields = [subfields;gridfields];
index =  eventdata.Indices(1);

% -- check which field was selected
if index > size(fields,1)-24 % last 24 are Gridjobparameters
    if ~ischar(allJobs{selected}.params.Gridjob.(fields{eventdata.Indices(1)}))
        elem = eval(eventdata.EditData); %evaluate string data to actual types
    else
        elem = eventdata.NewData; %evaluate strings as strings
    end
    
    allJobs{selected}.params.Gridjob.(fields{eventdata.Indices(1)}) = elem; %change parameter in current job-object
    %-- if jobname is changed, update listbox
    if index == size(fields,1)-15
        current = cellstr(get(handles.listbox3,'String'));
        current{selected+1} = allJobs{selected}.params.Gridjob.jobname;
        set(handles.listbox3,'String',current);
    end
else
    if ~ischar(allJobs{selected}.params.(val).(fields{eventdata.Indices(1)}))
        elem = eval(eventdata.EditData);
    else
        elem = eventdata.NewData;
    end
    allJobs{selected}.params.(val).(fields{eventdata.Indices(1)}) = elem;
end


% --- Executes on button press in 'add job' button.
% --- Allows to add job for execution
function Add_job_Callback(hObject, eventdata, handles)
% hObject    handle to Add_job (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global allJobs;
global num_jobs;
global selected;
global constructedFromFolder;

list = get(handles.popupmenu1,'String'); %classes
val = list{get(handles.popupmenu1,'Value')};%get chosen class
make = str2func(val(1:end-2));
obj = make();
current = cellstr(get(handles.listbox3,'String'));
if ~ismember(current,obj.params.Gridjob.jobname) %don't add job with same name (maybe remove, to allow copy of job)

    set(handles.listbox3,'String',[current;{obj.params.Gridjob.jobname}]); %add new job to listbox

    num_jobs = num_jobs + 1; %increase number of jobs
    selected = num_jobs;
    allJobs{num_jobs} = obj; %store job-object

    if ~isempty(constructedFromFolder)
        allJobs{num_jobs}.params.Gridjob.relativeWorkpath = char(constructedFromFolder); %set workpath as selected
    else
        workdir_Callback(hObject,eventdata,handles); %if none selected open dialog
        allJobs{num_jobs}.params.Gridjob.relativeWorkpath = char(constructedFromFolder); %set workpath as selected
    end


    set(handles.table1,'ColumnEditable',[false true false]);
end



% --- Executes on selection change in listbox3.
% --- Allows to reopen already added job
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3
global allJobs;
global selected;
get(handles.figure1,'SelectionType');
if strcmp(get(handles.figure1,'SelectionType'),'open') %find double click
    
    selected = get(handles.listbox3,'Value');
    selected = selected -1;
    obj = allJobs{selected};
    
    %-- update parameter table
    name = fieldnames(obj.params);
    name = name{2};
    [fields, values] = build_table(hObject,handles,obj,name);
    set(handles.table1,'data',[fields, values']);
    
    %set popupmenu and make value column editable
    str = get(handles.popupmenu1,'String');
    set(handles.popupmenu1,'Value',find(ismember(str,strcat(name,'.m'))));
    set(handles.table1,'ColumnEditable',[false true false]);
end

% --- Executes during object creation, after setting all properties.
% --- Listbox of added jobs
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in the 'remove' button
% --- Allows to remove job from listbox
function Remove_Callback(hObject, eventdata, handles)
% hObject    handle to Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global allJobs;
global selected;
global num_jobs;

sel = get(handles.listbox3,'Value');
prev_str = get(handles.listbox3, 'String');
if sel ~= 1
    if ~isempty(prev_str)
        prev_str(get(handles.listbox3,'Value')) = [];
        set(handles.listbox3, 'String', prev_str, ...
            'Value', min(sel,length(prev_str)));
        allJobs(sel-1) = []; %remove job
        num_jobs = num_jobs - 1; %decrease number of jobs
        selected = 0; %deselect
        popupmenu1_Callback(hObject, eventdata, handles);
    end
end

% --- Executes on button press in Add_param.
% --- Adds Parameter to object
function Add_param_Callback(hObject, eventdata, handles)
% hObject    handle to Add_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global allJobs;
global selected;

editable = get(handles.table1,'Columneditable'); %only if job can be edited / is selected
if editable(2)
    parameter = get(handles.edit1,'String');
    value = get(handles.edit2,'String');
    if ~strcmp('New parameter',parameter) && ~strcmp('New value',value)
        
        %set parameter
        param_parts = strsplit(parameter,'.');
        
        result = 'allJobs{selected}.params.';
        for k=1:size(param_parts,2)
            result = strcat(result,'(''',char(param_parts(k)),''').');
        end
        result = strcat(result(1:end-1), ' = ', value,';');
        eval(result);
        % reload table
        list = get(handles.popupmenu1,'String');
        val = list{get(handles.popupmenu1,'Value')};
        val = val(1:end-2);
        [fields values] = build_table(hObject,handles,allJobs{selected},val);
        set(handles.table1,'data',[fields, values']);
    end
end


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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit1.
function edit1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Toggel the "Enable" state to ON

set(hObject, 'Enable', 'On');

% Create UI control

uicontrol(handles.edit1);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit2.
function edit2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Toggel the "Enable" state to ON

set(hObject, 'Enable', 'On');

% Create UI control

uicontrol(handles.edit2);


% --- Executes on button press in workdir.
% --- Allows to select relative working directory
function workdir_Callback(hObject, eventdata, handles)
% hObject    handle to workdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global constructedFromFolder;
paths = dataPaths( );

matching = false;
while ~matching
    selected_dir = uigetdir(paths.workdir);
    if selected_dir == 0
        response = warndlg('Path has to be working directory!');
        uiwait(response);
        continue;
    end
    % save folder from which the object is constructed
    constructedFromFolder = pwd;
    paths = dataPaths( );
    %find relative path from pwd to paths.paramdir
    filesepTmp=filesep();
    if filesepTmp=='\'
        filesepTmp='\\';
    end
    partsParamdir = strread(fullfile(paths.workdir),'%s','delimiter',filesepTmp);
    partsPwd = strread(selected_dir,'%s','delimiter',filesepTmp);
    matching = true;
    for k=1:min(length(partsParamdir),length(partsPwd))
        if ~strcmp(partsParamdir{k},partsPwd{k})
            matching = false;
            break;
        end
    end
    if matching
        constructedFromFolder = partsPwd(length(partsParamdir)+1:end);
        constructedFromFolder = fullfile(constructedFromFolder{:});
    else
        constructedFromFolder = selected_dir;
        response = warndlg('Path has to be working directory!');
        uiwait(response);
    end
end


% --- Executes on button press in load.
% --- Allows to load job from file
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global allJobs;
global num_jobs;
global selected;
num_jobs = num_jobs+1;
selected = num_jobs;

paths = dataPaths( );
[file,selected_job] = uigetfile(paths.workdir); %get file
allJobs{num_jobs} = Gridjob(strcat(selected_job,file)); %create job object

names = fieldnames(allJobs{selected}.params);
val = names{2}; %get class type of object (names{1} = Gridjob)
[fields, values] = build_table(hObject,handles,allJobs{selected},val);
set(handles.table1,'data',[fields, values'])

str = get(handles.popupmenu1,'String');
set(handles.popupmenu1,'Value',find(ismember(str,strcat(val,'.m'))));
current = cellstr(get(handles.listbox3,'String'));
set(handles.listbox3,'String',[current;{allJobs{num_jobs}.params.Gridjob.jobname}]);
set(handles.table1,'ColumnEditable',[false true false]);


% --- Executes on button press in stop.
% --- Stop refresh of status table
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global running;
global allJobs;
running = ~running;
if running
    show_status(hObject,handles,allJobs);
end


% --- Executes when selected cell(s) is changed in status.
% --- Opens error- and logfiles on mouse click
function status_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to status (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
global constructedFromFolder;
c = config;
data = get(handles.status,'data');
if eventdata.Indices
    if eventdata.Indices(2) == 6 && ~isempty(data{eventdata.Indices(1),6})
        filename = fullfile(c.work_path,constructedFromFolder,data(eventdata.Indices(1),6));
        if exist(filename{1},'file')
            open(filename{1});
        end
    elseif eventdata.Indices(2) == 7 && ~isempty(data{eventdata.Indices(1),7})
        filename = fullfile(c.work_path,constructedFromFolder,data(eventdata.Indices(1),7));
        if exist(filename{1},'file')
            open(filename{1});
        end
    end
end


% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global allJobs;
global selected;
if ismethod(allJobs{selected},'plotJob') && finished(allJobs{selected}) %if job was loaded its of type Gridjob -> does not find method
    plotJob(allJobs{selected});
else
end

function finish = finished(job)
global constructedFromFolder;

con = config;
name = job.params.Gridjob.jobname;
filename = fullfile(con.work_path,constructedFromFolder,strcat('temp_',char(name)))
allFiles = dir(filename{1});
allNames = {allFiles.name};
%look for file isfinished - if it is deleted the job is finished
if ~ismember(allNames,'isfinished')
    finish = true;
else
    finish = false;
end

