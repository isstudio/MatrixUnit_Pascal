uses matrixunit_dev;

begin
  var m:TMatrix;
  m:= TMatrix.Create(3, 3, true);
  m.FillByKeyboard(false);
  m.ShowUp(false);
end.