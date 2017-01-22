{

    |\      /|  /-------\  |-----\    |-------|
    | \    / |  /       \  |      \   |
    |  \  /  |  /       \  |       \  |
    |   \/   |  /-------\  |       |  |----|
    |        |  /       \  |       /  |
    |        |  /       \  |      /   |
    |        |  /       \  |-----/    |-------|

                ---------  | /----\
                    |      |/      \
                    |      |       |
                    |      |       |
                    |      |       |
                    |      |       |
                ---------  |       |

     ---------  /-------\  /-------\  \-------/
         |      |          |           \-----/
         |      |          |            \---/
         |      \-------\  \-------\     \-/
         |              |          |      |
         |              |          |     |-|
     ---------  \-------/  \-------/     |-|
     
     /-------\  /-------\     /|      /-------\
     |       |  |       |    / |      |      
             |  |       |   /  |      |     
     /-------/  |       |  /   |      /-------\
     |          |       |      |      |       |   
     |          |       |      |      |       |   
     \-------/  \-------/  \-------/  \-------/  
      
}

(* WARNING!!! It's OLD version
   we don't use it now!!!      *)

unit DetermineUnit;

interface



const error_msg = 'Ошибка при работе с матрицей!';
      error_det = 'Детерминант данной матриц невозможно вычислить';

type
  int = longint;
  float = real;
  table = aRRay of aRRay of float;
  list = aRRay of float;
  TMatrix = class
  private
    //mas:table;
    function Minor(a:table; x, y:Byte):table;
    function pow(a, n:int):int;
  public
    mas:table;
    //value:list;
    constructor Create(x, y:Byte);//+
    procedure SetLen(x, y:Byte);//+
    procedure FillByKeyboard;//+
    procedure FillByFile(var fin:TextFile);//+
    procedure FillByMatrix(fin:table);//+
    procedure ShowUp;//+
    procedure SaveToFile(var fout:TextFile);//+
    procedure SaveToMatrix(var fout:table);//+
    function Determinant(a:table):float;
    //function KramerSolve:list;
    //function GaussSolve:list;
  end;
  
implementation

const k = -1;

constructor TMatrix.Create(x, y:Byte);
var i, j:Byte;
begin
  inherited Create;
  SetLength(mas, x);
  for i:= low(mas) to high(mas) do
  begin
    SetLength(mas[i], y);
    for j:= low(mas[i]) to high(mas[i]) do
      mas[i, j]:= 0;
  end;
end;

function TMatrix.pow(a, n:int):int;
begin
  if (n = 1) then
    Result:= a
  else Result:= a * pow(a, n - 1);
end;

procedure TMatrix.SetLen(x, y:Byte);
var i:Byte;
begin  
  SetLength(mas, x);
  for i:= low(mas) to high(mas) do
    SetLength(mas[i], y);
end;

procedure TMatrix.FillByKeyboard;
var i, j:Byte;
begin
  for i:= low(mas) to high(mas) do
    for j:= low(mas[i]) to high(mas[i]) do
      try
        read(mas[i, j]);
      except
        writeln(error_msg);
        exit;
      end;
end;

procedure TMatrix.FillByFile(var fin:TextFile);
var i, j:Byte;
begin
  try
    reset(fin);
    for i:= low(mas) to high(mas) do
      for j:= low(mas[i]) to high(mas[i]) do
        if eof(fin) then
        begin
          writeln(error_msg);
          exit;
        end
        else
        try
          read(fin, mas[i, j]);
        except
          writeln(error_msg);
          exit;
        end;
  finally
    closefile(fin);
  end;
end;

procedure TMatrix.FillByMatrix(fin:table);
var i, j:Byte;
begin
  for i:= low(mas) to high(mas) do
    for j:= low(mas[i]) to high(mas[i]) do
      try
        mas[i, j]:= fin[i, j];
      except
        writeln(error_msg);
        exit;
      end;
end;

procedure TMatrix.ShowUp;
var i, j:Byte;
begin
  for i:= low(mas) to high(mas) do
  begin
    for j:= low(mas[i]) to high(mas[i]) do
      write(mas[i, j]:4);
    writeln;
  end;
end;

procedure TMatrix.SaveToFile(var fout:TextFile);
var i, j:Byte;
begin
  try
    rewrite(fout);
    for i:= low(mas) to high(mas) do
    begin
      for j:= low(mas[i]) to high(mas[i]) do
        write(fout, mas[i, j]:4);
      writeln(fout);
    end;
  finally
    closefile(fout);
  end;
end;

procedure TMatrix.SaveToMatrix(var fout:table);
var i, j:Byte;
begin
  SetLength(fout, Length(mas));
  for i:= low(fout) to high(fout) do
  begin
    SetLength(fout[i], Length(mas[i]));
    for j:= low(fout[i]) to high(fout[i]) do
      fout[i, j]:= mas[i, j];
  end;
end;

function TMatrix.Determinant(a:table):float;
var det, max, tmp:float;
    i, j, n, y, u:Byte;
    buf:list;
    b:table;
begin
  i:= Length(a);
  j:= Length(a[low(a)]);
  if (i <> j) then
  begin
    writeln('Невозможно вычислить определитель');
    exit;
  end;
  if (i = 2) then
  begin
    Result:= a[0, 0] * a[1, 1] - a[1, 0] * a[0, 1];
    exit;
  end;
  det:= 1;
  n:= low(a);
  max:= abs(a[n, low(a[n])]);
  for i:= low(a) to high(a) do
  if (abs(a[i, low(a[i])]) > max) then
  begin
    max:= a[i, low(a[i])];
    n:= i;
  end;
  if (max = 0) then
  begin
    Result:= 0;
    exit;
  end;
  if (n <> low(a)) then
  begin
    buf:= a[n];
    a[n]:= a[low(a)];
    a[low(a)]:= buf;
    det:= -det;
  end;
  max:= a[low(a), low(a[low(a)])];
  for j:= low(a[low(a)]) to high(a[low(a)]) do
  begin
    a[low(a), j]:= a[low(a), j] / max;
    if (j > low(a[low(a)])) then
    begin
      tmp:= a[low(a), j];
      for i:= low(a) to high(a) do
        a[i, j]:= a[i, j] - tmp * a[i, low(a[i])];
    end;
  end;
  det:= det * max;
  SetLength(b, Length(a) - 1);
  for i:= low(b) to high(b) do
  begin
    SetLength(b[i], Length(a[i]) - 1);
    for j:= low(b[i]) to high(b[i]) do
    begin
      b[i, j]:= a[i + 1, j + 1];
    end;
  end;
  Result:= det * Determinant(b);
end;
  
{var i, j:Byte;
    rez:int;
begin
  i:= Length(a);
  j:= Length(a[0]);
  if not(i = j) then
  begin
    writeln(error_det);
    exit;
  end;
  if (i = 2) then
  begin
    writeln('found minimal matrix');
    Result:= a[0, 0] * a[1, 1] - a[1, 0] * a[0, 1]
  end
  else
  begin
    rez:= 0;
    for j:= low(a[0]) to high(a[0]) do
    begin
      writeln(a[0, j], pow(k, j));  
      inc(rez, a[0, j] * pow(k, j) * Determinant(Minor(a, 0, j)));
    end;
    Result:= rez;
  end;
end;}

function TMatrix.Minor(a:table; x, y:Byte):table;
var i, j, m:Byte;
    b:table;
    dx, dy:0..1;
begin
  i:= Length(a);
  m:= Length(a[0]);
  SetLength(b, i - 1);
  for i:= low(a) to high(a) do
  begin
    SetLength(mas[i], m - 1);
    for j:= low(a[i]) to high(a[i]) do
    begin
      if (i >= x) then
        dx:= 1
      else
        dx:= 0;
      if (j >= y) then
        dy:= 1
      else
        dy:= 0;
      b[i, j]:= a[i + dx, j + dy];
    end;
  end;
  for i:= low(b) to high(b) do
  begin
    for j:= low(b[i]) to high(b[i]) do
      write(b[i, j]:4);
    writeln;
  end;
  Result:= b;
end;
end.