---
title: FPGA常用基础模块
tags: [FPGA,Verilog]
categories: FPGA
index_img: /img/FPGA.jpg
banner_img: /img/default.png
---



数字电路中常用基础模块

<!-- more -->

# FPGA常用基础模块

列举一些比较典型的数字电路，包括组合逻辑电路和时序逻辑电路。体现一定的Verilog语法与编程技巧。同时这些模块也是作为实际电路设计中比较通用的模块。

## 组合逻辑电路

### 1.  8-1数据选择器

```verilog
module choose(
    input [2:0] A, 
    input [7:0] D,
    output Y
    );
reg y_temp;
always @(A or D)
begin
		case(A)
			3'b000: y_temp = D[0];
			3'b001: y_temp = D[1];
			3'b010: y_temp = D[2];
			3'b011: y_temp = D[3];
			3'b100: y_temp = D[4];
			3'b101: y_temp = D[5];
			3'b110: y_temp = D[6];
			default: y_temp = D[7];
		endcase
end
assign Y = y_temp;  
endmodule
```

### 2. 3-8译码器

```verilog
module Tran_3_8(
    input [2:0] In,
    input EN, //使能信号
    output [7:0] Out
    );
reg [7:0] Out_temp;
always @(In or EN)
begin
	if(~EN)	
		case(In)
			3'b000: Out_temp = 8'b0111_1111;
			3'b001: Out_temp = 8'b1011_1111;
			3'b010: Out_temp = 8'b1101_1111;
			3'b011: Out_temp = 8'b1110_1111;
			3'b100: Out_temp = 8'b1111_0111;
			3'b101: Out_temp = 8'b1111_1011;
			3'b110: Out_temp = 8'b1111_1101;
			3'b111: Out_temp = 8'b1111_1110;
			default: Out_temp = 8'b1111_1111;
		endcase
	else
		Out_temp = 8'b1111_1111;	
end                                                               
assign Out = Out_temp;
endmodule
```



## 时序逻辑电路

### 1. D触发器

```verilog
module DFF(S,R,D,CLK,Q,qn);
    
input S,R,D,CLK;
output reg Q;
output reg qn;
	
always @(posedge CLK or negedge S or negedge R)
begin
	if(!R) begin Q <= 1'b0; qn <= 1'b1; end
	else if(!S) begin Q <= 1'b1; qn <= 1'b0; end
	else begin Q <= D; qn <= ~D; end
end
endmodule
```

### 2. 模M同步二进制加法计数器

```verilog
module cnt_M(clk,rst,cnt);

parameter M = 1024;        //M mo
parameter N = 10;          //N wei count

input clk,rst;
output reg [N-1:0]cnt;

//M wei 10 de count
always @(posedge clk or negedge rst)
begin
	if(!rst) cnt <= 1'b0;
	else if(cnt < M-1) cnt <= cnt + 1'b1;
	else cnt <= 1'b0;
end

endmodule
```

### 3. 分频电路

```verilog
module Div_Clk(clk,rst,led);

input clk, rst;
output Div; 

reg  [25:0]  cnt; //计数值
parameter CNT = 32'd10;  //周期为CNT
parameter Half_CNT = 32'd5;  //周期的一半
//实际为一个计数器    
always @(posedge clk or negedge rst)
begin
    if (!rst) cnt <= 0;
    else if (cnt == CNT-1)  cnt<= 0;
    else  cnt <= cnt + 1'b1; 
end
    
assign led = (cnt < Half_CNT-1) ? 1'b1 : 1'b0;
    
endmodule
```

### 4. 有限状态机

状态机一般包括组合逻辑电路和寄存器两部分。

状态机的下一个状态的输出不仅与输入信号有关，还与寄存器当前状态有关。

状态机可以分为米勒（Mealy）型和摩尔（Moore）型。

Mealy型状态机的输出是当前状态和输入信号的函数，Moore型状态机的输出仅是当前状态的函数。

下面是一个三段式Moore型状态机。

在编写一个状态机之前，首先要绘制对应时序逻辑的状态转移图，然后根据状态转移图，编写程序。

```verilog
module TOP(rst,clk,yout);

input rst,clk;
output reg [1:0] yout;

parameter s0 = 3'b100, s1 = 3'b010, s2 = 3'b001;//定义状态机中的所有状态
reg [2:0] state;
reg [2:0] next_state;

//第一个always块，完成状态转换
always @(posedge clk or negedge rst)
begin
	if(!rst)  state <= s0;
	else state <= next_state;
end
//第二个always块，完成状态机的内部逻辑
always @(state or next state)
begin
    case(state)
        s0: next_state <= s1;
        s1: next_state <= s2;
        s2: next_state <= s0;
        default: next_state <= s0;
    endcase
end
//状态机的外部输出
always @(*)
begin
    case(state)
        s0: yout <= 2'b00;
        s1: yout <= 2'b01;
        s2: yout <= 2'b10;
        default: yout <= 2'b00;
    endcase
end

endmodule
```



