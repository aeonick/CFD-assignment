classdef Solver<handle
%- Author:石凯元
%- Time: 08 Jul 2019
%- Please follow GPL License using the source code
    properties
        mesh;model;enableGPU;
    end
    methods
        function self = Solver(mesh,model,enableGPU)
            self.mesh=mesh;
            self.model=model;
            self.enableGPU=enableGPU;
        end
        function [array]=gpuArray(self,array)
            if self.enableGPU
                array=gpuArray(array);
            end
        end
        function gpuIn(self)
            if self.enableGPU==0
                return
            end
            self.model.u=self.gpuArray(self.model.u);
            self.model.v=self.gpuArray(self.model.v);
            self.model.p=self.gpuArray(self.model.p);
            names=fieldnames(self.mesh.operators);
            for i=1:length(names)
                name=names{i};
                operator=getfield(self.mesh.operators,name);
                if strcmp(class(operator),'Operator')
                    operator.gpuIn();
                    self.mesh.operators=setfield(self.mesh.operators,name,operator);
                else
                    self.mesh.operators=setfield(self.mesh.operators,name,self.gpuArray(operator));
                end
            end
        end
        function gpuOut(self)
            self.model.u=gather(self.model.u);
            self.model.v=gather(self.model.v);
            self.model.p=gather(self.model.p);
            names=fieldnames(self.mesh.operators);
            for i=1:length(names)
                name=names{i};
                operator=getfield(self.mesh.operators,name);
                if strcmp(class(operator),'Operator')
                    operator.gpuOut();
                    self.mesh.operators=setfield(self.mesh.operators,name,operator);
                else
                    self.mesh.operators=setfield(self.mesh.operators,name,gather(operator));
                end
            end
        end
        function NSSolver(self,dt,steps,ifOptimize)
            self.mesh.createOperator();
            self.gpuIn();
            mesh=self.mesh;
            model=self.model;
            u=model.u;
            v=model.v;
            p=model.p;
            b = self.gpuArray(zeros(mesh.num,1));
            if ifOptimize
                mesh.optimize();
            end
            op=mesh.operators;
            for n=1:steps
                un = u;
                vn = v;
                b=(1/dt*(op.CDX*u+op.CDY*v)*(mesh.dx*mesh.dy)-(op.CDX*u).^2-2*(op.CDY*u).*(op.CDX*v)-2*(op.CDY*v).*(op.CDX*u))*model.rho;
                for q=1:4
                    p=(op.int*p-b)/(2*(mesh.dx^2+mesh.dy^2));
                    p=p+op.neumann*p;
                end
                temp=model.nu.*op.laplace-mesh.UpWindX(un)-mesh.UpWindY(vn);
                u=un+(temp*un-1/model.rho*(op.CDX*p))*(dt/mesh.dx/mesh.dy);
                v=vn+(temp*vn-1/model.rho*(op.CDY*p))*(dt/mesh.dx/mesh.dy);
            end
            self.model.u=u;
            self.model.v=v;
            self.model.p=p;
            self.gpuOut();
        end
        function dataMesh=show(self,data)
            dataMesh=self.mesh.show(data);
        end
    end
end