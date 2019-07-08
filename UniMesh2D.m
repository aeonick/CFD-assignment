classdef UniMesh2D<handle
%- Author:石凯元
%- Time: 08 Jul 2019
%- Please follow GPL License using the source code
    properties
        x;y;dx;dy;
        row;col;map;
        num;toler;
        e;w;n;s;c;
        indexs;fields;
        operators;
        boundrys;
    end
    methods
        function self = UniMesh2D(varargin)
            if isempty(varargin)
                varargin={[],[]};
            end
            self.x=varargin{1}(:);
            self.y=varargin{2}(:);
            self.num=length(self.x);
            self.toler=1e-14;
        end
        function newObj = plus(self1,self2)
            x=[self1.x;self2.x];
            y=[self1.y;self2.y];
            uni=uniquetol([x y],1e-14,'ByRows',true);
            newObj = UniMesh2D(uni(:,1),uni(:,2));
            newObj.createIndex();
        end
        function createBlock(self,x_start,x_end,y_start,y_end,num_x,num_y)
            x=linspace(x_start,x_end,num_x);
            y=linspace(y_start,y_end,num_y);
            [x,y]=meshgrid(x,y);
            x=x';y=y';
            self.x=x(:);self.y=y(:);
            self.num=length(self.x);
            self.createIndex();
        end
        function createIndex(self)
            self.c=(1:self.num)';
            self.e=zeros(self.num,1);
            self.w=zeros(self.num,1);
            self.n=zeros(self.num,1);
            self.s=zeros(self.num,1);
            rowMap=unique(self.x,'sorted');
            colMap=unique(self.y,'sorted');
            self.dx=rowMap(2)-rowMap(1);
            self.dy=colMap(2)-colMap(1);
            for i=1:length(rowMap)
                self.row(abs(self.x-rowMap(i))<self.toler)=i;
            end
            for i=1:length(colMap)
                self.col(abs(self.y-colMap(i))<self.toler)=i;
            end
            self.map=zeros(length(rowMap)+2,length(colMap)+2);
            self.map=self.addData(self.map,self.row+1,self.col+1,1:self.num);
            for i=1:self.num
                self.e(i)=findDir(self,i,1,0);
                self.w(i)=findDir(self,i,-1,0);
                self.n(i)=findDir(self,i,0,1);
                self.s(i)=findDir(self,i,0,-1);
            end
        end
        function data=addData(self,data,row,col,value)
            data(sub2ind(size(data),row,col))=data(sub2ind(size(data),row,col))+value;
        end
        function re=findDir(self,i,dirX,dirY)
            r=self.row(i)+dirX+1;
            c=self.col(i)+dirY+1;
            if self.map(r,c)==0
                r=self.row(i)-dirX+1;
                c=self.col(i)-dirY+1;
            end
            re=self.map(r,c);
        end
        function dataMesh=show(self,t)
            dataMesh=double(self.map);
            dataMesh=dataMesh(2:end-1,2:end-1);
            index=find(dataMesh);
            dataMesh(dataMesh==0)=NaN;
            dataMesh(index)=t(dataMesh(index));
        end
        function createOperator(self)
            dx=self.dx;dy=self.dy;
            laplace=zeros(self.num);
            laplace=self.addData(laplace,self.c,self.c,-2*dx/dy);
            laplace=self.addData(laplace,self.c,self.n,1*dx/dy);
            laplace=self.addData(laplace,self.c,self.s,1*dx/dy);
            laplace=self.addData(laplace,self.c,self.c,-2*dy/dx);
            laplace=self.addData(laplace,self.c,self.w,1*dy/dx);
            laplace=self.addData(laplace,self.c,self.e,1*dy/dx);
            laplace(self.boundrys.dirichlet,:)=0;
            self.operators.laplace=laplace;
            centerDiffY=zeros(self.num);
            centerDiffY=self.addData(centerDiffY,self.c,self.n,0.5*dx);
            centerDiffY=self.addData(centerDiffY,self.c,self.s,-0.5*dx);
            centerDiffY(self.boundrys.dirichlet,:)=0;
            self.operators.CDY=centerDiffY;
            centerDiffX=zeros(self.num);
            centerDiffX=self.addData(centerDiffX,self.c,self.e,0.5*dy);
            centerDiffX=self.addData(centerDiffX,self.c,self.w,-0.5*dy);
            centerDiffX(self.boundrys.dirichlet,:)=0;
            self.operators.CDX=centerDiffX;
            unlinDiffX0=zeros(self.num);
            unlinDiffX0=self.addData(unlinDiffX0,self.c,self.w,-0.5*dy);
            unlinDiffX0=self.addData(unlinDiffX0,self.c,self.c,1*dy);
            unlinDiffX0=self.addData(unlinDiffX0,self.c,self.e,-0.5*dy);
            unlinDiffX0(self.boundrys.dirichlet,:)=0;
            self.operators.UDX0=unlinDiffX0;
            unlinDiffX1=zeros(self.num);
            unlinDiffX1=self.addData(unlinDiffX1,self.c,self.e,0.5*dy);
            unlinDiffX1=self.addData(unlinDiffX1,self.c,self.w,-0.5*dy);
            unlinDiffX1(self.boundrys.dirichlet,:)=0;
            self.operators.UDX1=unlinDiffX1;
            unlinDiffY0=zeros(self.num);
            unlinDiffY0=self.addData(unlinDiffY0,self.c,self.n,-0.5*dx);
            unlinDiffY0=self.addData(unlinDiffY0,self.c,self.c,1*dx);
            unlinDiffY0=self.addData(unlinDiffY0,self.c,self.s,-0.5*dx);
            unlinDiffY0(self.boundrys.dirichlet,:)=0;
            self.operators.UDY0=unlinDiffY0;
            unlinDiffY1=zeros(self.num);
            unlinDiffY1=self.addData(unlinDiffY1,self.c,self.n,0.5*dx);
            unlinDiffY1=self.addData(unlinDiffY1,self.c,self.s,-0.5*dx);
            unlinDiffY1(self.boundrys.dirichlet,:)=0;
            self.operators.UDY1=unlinDiffY1;
            int=zeros(self.num);
            int=self.addData(int,self.c,self.n,dx^2);
            int=self.addData(int,self.c,self.s,dx^2);
            int=self.addData(int,self.c,self.w,dy^2);
            int=self.addData(int,self.c,self.e,dy^2);
            int(self.boundrys.dirichlet,:)=0;
            self.operators.int=int;
        end
        function operator=UpWindX(self,coeff)
            operator=abs(coeff).*self.operators.UDX0+coeff.*self.operators.UDX1;
        end
        function operator=UpWindY(self,coeff)
            operator=abs(coeff).*self.operators.UDY0+coeff.*self.operators.UDY1;
        end
        function optimize(self)
            names=fieldnames(self.operators);
            for i=1:length(names)
                name=names{i};
                self.operators=setfield(self.operators,name,Operator(self.c,{self.c,self.n,self.e,self.w,self.s},getfield(self.operators,name)));
            end
        end
    end
end