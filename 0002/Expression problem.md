---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Julia 1.7.0-DEV
    language: julia
    name: julia-1.7
---

<!-- #region -->
# Expression problem

* https://miguelraz.github.io/blog/dispatch/index.html

以下のコードは

* https://eli.thegreenplace.net/2016/the-expression-problem-and-its-solutions/

より。

__C++__

```C++
class Expr {
public:
  virtual std::string ToString() const = 0;
  virtual double Eval() const = 0;
};
```

```C++
class Constant : public Expr {
public:
  Constant(double value) : value_(value) {}

  std::string ToString() const {
    std::ostringstream ss;
    ss << value_;
    return ss.str();
  }

  double Eval() const {
    return value_;
  }

private:
  double value_;
};
```

```C++
class BinaryPlus : public Expr {
public:
  BinaryPlus(const Expr& lhs, const Expr& rhs) : lhs_(lhs), rhs_(rhs) {}

  std::string ToString() const {
    return lhs_.ToString() + " + " + rhs_.ToString();
  }

  double Eval() const {
    return lhs_.Eval() + rhs_.Eval();
  }

private:
  const Expr& lhs_;
  const Expr& rhs_;
};
```

__Haskell__

```Haskell
module Expressions where

data Expr = Constant Double
          | BinaryPlus Expr Expr

stringify :: Expr -> String
stringify (Constant c) = show c
stringify (BinaryPlus lhs rhs) = stringify lhs
                                ++ " + "
                                ++ stringify rhs

evaluate :: Expr -> Double
evaluate (Constant c) = c
evaluate (BinaryPlus lhs rhs) = evaluate lhs + evaluate rhs
```
<!-- #endregion -->

## 値の四則演算を行えるモジュールO

* 式の型
  * 定数の型
  * 二項演算達の型
* これらを文字列に変換する函数
* そららの値を評価する函数

```julia
module O

abstract type Expression end

struct Constant{C} <: Expression value::C end
value(x::Constant) = getfield(x, :value)
stringify(x::Constant) = string(value(x))
evaluate(x::Constant) = value(x)

abstract type BinOp <: Expression end
lhs(x::BinOp) = getfield(x, :lhs)
rhs(x::BinOp) = getfield(x, :rhs)

const binop_list = ((:Plus, :+), (:Minus, :-), (:Mult, :*), (:Div, :/))
for (S, op) in binop_list
    @eval begin
        struct $S{L<:Expression, R<:Expression} <: BinOp lhs::L; rhs::R end
        stringify(x::$S) = "(" * stringify(lhs(x)) * " $($op) " * stringify(rhs(x)) * ")"
        evaluate(x::$S) = $op(evaluate(lhs(x)), evaluate(rhs(x)))
    end
end

# デフォルトでの表示の仕方
Base.show(io::IO, x::Expression) = print(io, stringify(x))

end
```

```julia
a = O.Constant(2)
b = O.Constant(3)
c = O.Constant(4)
d = O.Constant(5)
e = O.Constant(6)
a, b, c, d, e
```

```julia
expr1 = O.Div(O.Plus(O.Minus(c, a), O.Mult(b, e)), d)
```

```julia
O.evaluate(expr1)
```

## モジュールOを変更することなく、モジュールPで型と函数を追加

モジュールPでは以下をO.Expression型に追加する。

* 函数呼び出しの型
* 函数呼び出しを文字列に変換するメソッド
* 函数を呼び出すメソッド
* 変数の型
* 変数を文字列に変換するメソッド
* 変数の値のリストの型
* 変数の値のリストを文字列に変換するメソッド
* 変数の値のリストに従って式を評価する函数
* 変数の値のリストに従って変数に式を代入した式を作る函数

```julia
module P

using InteractiveUtils: subtypes
using ..O: O, Expression, Constant, BinOp, lhs, rhs

struct FunctionCall{F, X<:Expression} <: Expression f::F; x::X end
fun(fx::FunctionCall) = getfield(fx, :f)
arg(fx::FunctionCall) = getfield(fx, :x)
O.stringify(fx::FunctionCall) = string(fun(fx)) * "(" * O.stringify(arg(fx)) * ")"
O.evaluate(fx::FunctionCall) = fun(fx)(O.evaluate(arg(fx)))

struct Variable <:Expression name::Symbol end
name(x::Variable) = getfield(x, :name)
O.stringify(x::Variable) = string(name(x))

struct ValueList{T<:NamedTuple} v::T end
ValueList() = ValueList((;))
ValueList(; p...) = ValueList((; p...))
parent(x::ValueList) = getfield(x, :v)
names(x::ValueList) = keys(parent(x))
value(x::ValueList, name::Symbol) = getproperty(parent(x), name)
O.stringify(x::ValueList) = "ValueList" * string(parent(x))

evaluate(x::Variable, v::ValueList) = value(v, name(x))
evaluate(x::FunctionCall, v::ValueList) = fun(x)(evaluate(arg(x), v))
evaluate(x::Constant, v::ValueList) = O.value(x)
for (S, op) in O.binop_list
    @eval evaluate(x::O.$S, v::ValueList=ValueList()) = $op(evaluate(lhs(x), v), evaluate(rhs(x), v))
end

substitute(x::Variable, v::ValueList) = name(x) ∈ names(v) ? Constant(value(v, name(x))) : x
substitute(x::FunctionCall, v::ValueList) = FunctionCall(fun(x), substitute(arg(x), v))
substitute(x::Constant, v::ValueList) = x
for (S, op) in O.binop_list
    @eval substitute(x::O.$S, v::ValueList) = O.$S(substitute(lhs(x), v), substitute(rhs(x), v))
end

# デフォルトでの表示の仕方
Base.show(io::IO, x::ValueList) = print(io, O.stringify(x))

end
```

```julia
P.ValueList(u=2, v=3, w=4, x=5, y=6)
```

```julia
u = P.Variable(:u)
v = P.Variable(:v)
w = P.Variable(:w)
x = P.Variable(:x)
y = P.Variable(:y)
u, v, w, x, y
```

```julia
sinpi_x_div_y = P.FunctionCall(sinpi, O.Div(x, y))
```

```julia tags=[]
P.evaluate(sinpi_x_div_y, P.ValueList(x = 1, y = 6))
```

```julia tags=[]
P.substitute(sinpi_x_div_y, P.ValueList(x = 1, y = 6))
```

```julia
P.evaluate(expr1, P.ValueList())
```

```julia
expr2 = O.Div(O.Plus(O.Minus(w, u), O.Mult(v, y)), x)
```

```julia
P.evaluate(expr2, P.ValueList(u=2, v=3, w=4, x=5, y=6))
```

```julia
expr3 = P.substitute(expr2, P.ValueList(u=2, v=3, w=4))
```

```julia
expr4 = P.substitute(expr3, P.ValueList(x = 5, y = 6))
```

```julia
O.evaluate(expr4)
```

```julia
P.evaluate(expr3, P.ValueList(x = 5, y = 6))
```

```julia

```
