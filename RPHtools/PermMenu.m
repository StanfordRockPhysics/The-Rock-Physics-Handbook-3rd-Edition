function Menu

Class.name(1)={'(Porosity, Grain Size, Shape Factor) Inputs Type Models'};
Kclass.name={'BernabeE';'FredrichE';'KozCarmE';'ModKozCarm';'PandaLakeKCE';
        'RevilE';'WylGregE'};

Class.name(2)={'(Porosity, Sw-irreducible) Inputs Type Models'};
Tclass.name={'CoatDum';'Coates';'Timur';'Tixier'};

Class.name(3)={'Statistical Based Models'};
Sclass.name={'Bloch';'Owolabi';'PandaLake'};

% get the class of model
classnames=getfield(Class,'name');
[class,action] = listdlg(...
	'liststring',	classnames,...
	'promptstring',	'SELECT THE TYPE OF MODEL',...
	'listsize',		[450 80],...
	'fus',		10,...
	'ffs',		10,...
	'SelectionMode',	'single');
if (action==0)
	return;
end

%get the name of the model
switch class
	case 1, modelnames = getfield(Kclass,'name');
	case 2, modelnames = getfield(Tclass,'name');
	case 3, modelnames = getfield(Sclass,'name');
end
[model,action] = listdlg(...
	'liststring',	modelnames,...
	'promptstring',	'SELECT THE MODEL',...
	'listsize',		[150 130],...
	'fus',		10,...
	'ffs',		10,...
	'SelectionMode',	'single');
if (action==0)
	return;
end

% run the model
switch class
	case 1, feval(Kclass.name{model});
	case 2, feval(Tclass.name{model});
	case 3, feval(Sclass.name{model});
end
