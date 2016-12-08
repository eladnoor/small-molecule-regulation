# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 16:40:16 2015

@author: noore
"""

from pythonds.basic.stack import Stack
import re
import numpy as np

class BoolParser(object):
    
    PREC = {'(': 1, ')': 1, 'and' : 2, 'or' : 2}
    
    @staticmethod
    def isBoolVariable(token):
        return re.findall('^b\d+$', token) != []
    
    def __init__(self, expr):
        # pad all paranthese with spaces to help the tokenizer
        expr = expr.replace('(', ' ( ')
        expr = expr.replace(')', ' ) ')
        self.prefix = self.infix_to_prefix(expr)
        self.postfix = self.infix_to_postfix(expr)

    def infix_to_prefix(self, infixexpr):
        def invert_paranth(s):
            if s == '(':
                return ')'
            elif s == ')':
                return '('
            else:
                return s
                
        tokenList = map(invert_paranth, infixexpr.split())
        tokenList.reverse()
        postfixList = self.infix_to_postfix(' '.join(tokenList))
        postfixList.reverse()
        return postfixList
    
    def infix_to_postfix(self, infixexpr):
        opStack = Stack()
        postfixList = []
        tokenList = infixexpr.split()
        
        for token in tokenList:
            if BoolParser.isBoolVariable(token):
                # if token is a boolean variable, just append to list
                postfixList.append(token)
            elif token == '(':
                # start a new nested expression in the stack
                opStack.push(token)
            elif token == ')':
                # end the nested expression by moving the operators
                # from the stack to the list (i.e. in reverse order)
                topToken = opStack.pop()
                while topToken != '(':
                    postfixList.append(topToken)
                    topToken = opStack.pop()
            else:
                while (not opStack.isEmpty()) and \
                   (BoolParser.PREC[opStack.peek()] >= BoolParser.PREC[token]):
                      postfixList.append(opStack.pop())
                opStack.push(token)
    
        while not opStack.isEmpty():
            postfixList.append(opStack.pop())
        
        return postfixList

    def __str__(self):
        return BoolParser.printPrefix(self.prefix)[0]

    @staticmethod
    def printPrefix(tokenList, pos=0):
        if BoolParser.isBoolVariable(tokenList[pos]): 
            return tokenList[pos], pos+1
        else: # this is an operator
            s1, pos1 = BoolParser.printPrefix(tokenList, pos+1)
            s2, pos2 = BoolParser.printPrefix(tokenList, pos1)
            return tokenList[pos] + '(' + s1 + ', ' + s2 + ')', pos2

    @staticmethod
    def calcPrefix(tokenList, value_dict, pos=0, defval=0):
        if tokenList == []:
            return np.nan, 0
        if BoolParser.isBoolVariable(tokenList[pos]): 
            return value_dict.get(tokenList[pos], defval), pos+1
        else: # this is an operator
            val1, pos1 = BoolParser.calcPrefix(tokenList, value_dict, pos+1, defval)
            val2, pos2 = BoolParser.calcPrefix(tokenList, value_dict, pos1, defval)

            if np.isnan(val1):
                val = val2
            elif np.isnan(val2):
                val = val1
            elif tokenList[pos] == 'and':
                val = min(val1, val2)
            elif tokenList[pos] == 'or':
                val = val1 + val2
            
            return val, pos2

    def evaluate(self, value_dict, defval=0):
        return BoolParser.calcPrefix(self.prefix, value_dict, defval=defval)[0]

#expr = "((b2 or (b5 and b4)) or (b1 and b4) or (b1 and b2)) and (b5 or b5)"
#a = BoolParser(expr)
#print expr
#print a.prefix
#print a
#
#x = a.evaluate({'b1': 10, 'b2': 20, 'b3': 30, 'b4': 40, 'b5': 50, 'b6': 60})
#print x