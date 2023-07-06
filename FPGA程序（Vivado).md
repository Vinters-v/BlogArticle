# FPGA程序（Vivado)

## DA

```verilog
module DA(clk,rst,DAdata1,DAdata2,daclk1,dawrt1,daclk2,dawrt2,AD_Data_in);

input clk,rst;
input [13:0] AD_Data_in;
output daclk1,dawrt1,daclk2,dawrt2;
output [13:0] DAdata1;
output [13:0] DAdata2;

wire [9:0] cnt;
reg [13:0] DA;

cnt_M cnt_M_10 (.clk(clk), .rst(rst),.cnt(cnt));

blk_mem_gen_0 DAout_sin (
  .clka(clk),    // input wire clka
  .ena(rst),      // input wire ena
  .addra(cnt),  // input wire [9 : 0] addra
  .douta(DAdata2)  // output wire [13 : 0] douta
);

always @(posedge clk or negedge rst)
begin
	if(!rst) DA <= 1'b0; 
	else 
		begin
		 DA[13:0] <= ~AD_Data_in[13:0];
		end
end

assign DAdata1 = DA;

assign daclk1 = clk;
assign dawrt1 = clk;
assign daclk2 = clk;
assign dawrt2 = clk;

endmodule




; Made by voiue
; WuHan,China
; width =14
; depth =16384
memory_initialization_radix=10;
memory_initialization_vector=	
```



## cnt_M

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



## AD

```verilog
module AD(clk,rst,AD_Data_1,AD_Data_2,out_1,out_2,AD_clk_1,AD_clk_2);

input clk,rst;
input [11:0] AD_Data_1,AD_Data_2;
output reg [13:0] out_1,out_2;
output AD_clk_1,AD_clk_2;

reg [11:0] CH_1,CH_2;

assign AD_clk_1 = clk;
assign AD_clk_2 = clk;

always @(posedge clk or negedge rst)
begin
	if(!rst) CH_1 <= 0;
	else
		begin
			CH_1[0] <= AD_Data_1[11];
			CH_1[1] <= AD_Data_1[10];
			CH_1[2] <= AD_Data_1[9];
			CH_1[3] <= AD_Data_1[8];
			CH_1[4] <= AD_Data_1[7];
			CH_1[5] <= AD_Data_1[6];
			CH_1[6] <= AD_Data_1[5];
			CH_1[7] <= AD_Data_1[4];
			CH_1[8] <= AD_Data_1[3];
			CH_1[9] <= AD_Data_1[2];
			CH_1[10] <= AD_Data_1[1];
			CH_1[11] <= AD_Data_1[0];
		end
end

always @(posedge clk or negedge rst)
begin
	if(!rst) CH_2 <= 0;
	else
		begin
			CH_2[0] <= AD_Data_2[11];
			CH_2[1] <= AD_Data_2[10];
			CH_2[2] <= AD_Data_2[9];
			CH_2[3] <= AD_Data_2[8];
			CH_2[4] <= AD_Data_2[7];
			CH_2[5] <= AD_Data_2[6];
			CH_2[6] <= AD_Data_2[5];
			CH_2[7] <= AD_Data_2[4];
			CH_2[8] <= AD_Data_2[3];
			CH_2[9] <= AD_Data_2[2];
			CH_2[10] <= AD_Data_2[1];
			CH_2[11] <= AD_Data_2[0];
		end
end

always @(posedge clk)
begin
	out_1[13:2] <= CH_1[11:0];
	out_1[1] <= CH_1[11] & 1'b1;
	out_1[0] <= CH_1[11] & 1'b1;
end

always @(posedge clk)
begin
	out_2[13:2] <= CH_2[11:0];
	out_2[1] <= CH_2[11] & 1'b1;
	out_2[0] <= CH_2[11] & 1'b1;
end


endmodule

```



## LED

```verilog
//利用计数器实现周期为1s的呼吸灯
module LED(clk,rst,led);

input clk, rst;
output led;

reg  [25:0]  cnt;
//parameter CNT = 32'd50_000_000;
//parameter Half_CNT = 32'd25_000_000;

parameter CNT = 32'd10;
parameter Half_CNT = 32'd5;


always @(posedge clk or negedge rst)
begin
    if (!rst) cnt <= 0;
    else if (cnt == CNT)  cnt<= 0;
    else  cnt <= cnt + 1'b1; 
end

assign led = (cnt < Half_CNT) ? 1'b1 : 1'b0;

endmodule
```



## TOP

```verilog
module TOP(

//sys
input clk,
input rst,
//AD
input [11:0] AD_Data_1,
output AD_clk_1,
//DA
output [13:0] DAdata1,
output daclk1,
output dawrt1,
output [13:0] DAdata2,
output daclk2,
output dawrt2,
//led
output led

    );
    
//DA
//wire [13:0] DAdata2;
//wire daclk2;
//wire dawrt2;

//AD
wire AD_clk_2;
wire [11:0] AD_Data_2;
//mid_AD
wire [13:0] AD1_14;
wire [13:0] AD2_14;    
    
LED led_4 (
    .clk(clk), //wire
    .rst(rst),  
    .led(led)   //output
    );    
    
DA DAC9767 (
    .clk(clk), 
    .rst(rst), 
    .DAdata1(DAdata1),  //output
    .DAdata2(DAdata2),  //output
    .daclk1(daclk1),  //output
    .dawrt1(dawrt1),  //output
    .daclk2(daclk2),  //output
    .dawrt2(dawrt2),  //output
    .AD_Data_in(AD1_14)
    );

AD ADC9226 (
    .clk(clk),
    .rst(rst), 
    .AD_Data_1(AD_Data_1), //input
    .AD_Data_2(AD_Data_2), //input
    .out_1(AD1_14),    //wire
    .out_2(AD2_14),    //wire
	.AD_clk_1(AD_clk_1),   //output
    .AD_clk_2(AD_clk_2)    //output
    );    
    
endmodule

```



## Delay

```matlab
module Delay(
input clk,
input rst,
input [13:0]DataIn,
output reg [13:0]DataOut

    );

always @(posedge clk or negedge rst)
begin
    if(!rst) DataOut <= 0;
    else DataOut <= DataIn;
end 
 
endmodule
```



## sim

```verilog
module sim_1(
    );
    
reg clk;
reg rst;

wire led;

wire [13:0] DataIn;
wire [27:0] yout;

TOP TOP_Init (
    .clk(clk), 
    .rst(rst),  
    .DataIn(DataIn),
    .yout(yout),
    .led(led)   //output
    );
 
integer dout_file;
reg valid;

initial 
begin
    clk = 0;
    rst = 0;
    #50;
    rst = 1;
    valid = 1;
    dout_file = $fopen("1.txt","w");
    if(dout_file == 0)  begin
        $display ("can not open the file!");    //创建文件失败，显示can not open the file!
        $stop;
       end

end

always #100 clk = ~clk;

always @(negedge clk)
begin
     if(rst&&valid)  begin //使能     
       $fwrite(dout_file,"%d\n",yout);    //使能信号有效，每来一个时钟，写入到所创建的文件中
         //$fwrite(dout_file,"%d\n",$signed(yout));    //使能信号有效，每来一个时钟，写入到所创建的文件中
     end
     else begin
        $fclose(dout_file); 
     end
end
endmodule
```



## constrains

```verilog
set_property IOSTANDARD LVCMOS33 [get_ports rst]
set_property IOSTANDARD LVCMOS33 [get_ports clk]
set_property IOSTANDARD LVCMOS33 [get_ports led]

set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[13]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[12]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[11]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[10]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[9]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[8]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[7]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[6]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[5]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[4]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[3]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[2]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[1]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata2[0]}]

set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[13]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[12]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[11]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[10]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[9]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[8]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[7]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[6]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[5]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[4]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[3]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[2]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[1]}]
set_property IOSTANDARD LVCMOS33 [get_ports {DAdata1[0]}]

set_property IOSTANDARD LVCMOS33 [get_ports daclk1]
set_property IOSTANDARD LVCMOS33 [get_ports daclk2]
set_property IOSTANDARD LVCMOS33 [get_ports dawrt1]
set_property IOSTANDARD LVCMOS33 [get_ports dawrt2]

set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[11]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[10]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[9]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[8]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[7]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[6]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[5]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[4]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[3]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[2]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[1]}]
set_property IOSTANDARD LVCMOS33 [get_ports {AD_Data_1[0]}]
set_property IOSTANDARD LVCMOS33 [get_ports AD_clk_1]

set_property PACKAGE_PIN D18 [get_ports clk]
set_property PACKAGE_PIN A19 [get_ports led]
set_property PACKAGE_PIN A20 [get_ports rst]

set_property PACKAGE_PIN U2 [get_ports {DAdata1[13]}]
set_property PACKAGE_PIN R1 [get_ports {DAdata1[11]}]
set_property PACKAGE_PIN U6 [get_ports {DAdata1[9]}]
set_property PACKAGE_PIN T4 [get_ports {DAdata1[7]}]
set_property PACKAGE_PIN T2 [get_ports {DAdata1[5]}]
set_property PACKAGE_PIN P6 [get_ports {DAdata1[3]}]
set_property PACKAGE_PIN M6 [get_ports {DAdata1[1]}]
set_property PACKAGE_PIN N1 [get_ports dawrt1]

set_property PACKAGE_PIN U1 [get_ports {DAdata1[12]}]
set_property PACKAGE_PIN P1 [get_ports {DAdata1[10]}]
set_property PACKAGE_PIN U5 [get_ports {DAdata1[8]}]
set_property PACKAGE_PIN T3 [get_ports {DAdata1[6]}]
set_property PACKAGE_PIN R2 [get_ports {DAdata1[4]}]
set_property PACKAGE_PIN P5 [get_ports {DAdata1[2]}]
set_property PACKAGE_PIN M5 [get_ports {DAdata1[0]}]
set_property PACKAGE_PIN M1 [get_ports daclk1]

set_property PACKAGE_PIN M2 [get_ports daclk2]
set_property PACKAGE_PIN K1 [get_ports {DAdata2[13]}]
set_property PACKAGE_PIN L3 [get_ports {DAdata2[11]}]
set_property PACKAGE_PIN K3 [get_ports {DAdata2[9]}]
set_property PACKAGE_PIN H2 [get_ports {DAdata2[7]}]
set_property PACKAGE_PIN H9 [get_ports {DAdata2[5]}]
set_property PACKAGE_PIN F3 [get_ports {DAdata2[3]}]
set_property PACKAGE_PIN G2 [get_ports {DAdata2[1]}]

set_property PACKAGE_PIN L2 [get_ports dawrt2]
set_property PACKAGE_PIN J1 [get_ports {DAdata2[12]}]
set_property PACKAGE_PIN K2 [get_ports {DAdata2[10]}]
set_property PACKAGE_PIN J3 [get_ports {DAdata2[8]}]
set_property PACKAGE_PIN H1 [get_ports {DAdata2[6]}]
set_property PACKAGE_PIN G9 [get_ports {DAdata2[4]}]
set_property PACKAGE_PIN E3 [get_ports {DAdata2[2]}]
set_property PACKAGE_PIN G1 [get_ports {DAdata2[0]}]

set_property PACKAGE_PIN T25 [get_ports {AD_Data_1[11]}]
set_property PACKAGE_PIN W26 [get_ports {AD_Data_1[9]}]
set_property PACKAGE_PIN Y26 [get_ports {AD_Data_1[7]}]
set_property PACKAGE_PIN W24 [get_ports {AD_Data_1[5]}]
set_property PACKAGE_PIN AA25 [get_ports {AD_Data_1[3]}]
set_property PACKAGE_PIN AB25 [get_ports {AD_Data_1[1]}]

set_property PACKAGE_PIN T24 [get_ports AD_clk_1]
set_property PACKAGE_PIN V26 [get_ports {AD_Data_1[10]}]
set_property PACKAGE_PIN W25 [get_ports {AD_Data_1[8]}]
set_property PACKAGE_PIN V24 [get_ports {AD_Data_1[6]}]
set_property PACKAGE_PIN Y25 [get_ports {AD_Data_1[4]}]
set_property PACKAGE_PIN AA24 [get_ports {AD_Data_1[2]}]
set_property PACKAGE_PIN AB26 [get_ports {AD_Data_1[0]}]

```

