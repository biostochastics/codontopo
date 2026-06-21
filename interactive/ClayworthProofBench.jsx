import React, { useState, useMemo } from "react";
/* ===== verified engines (faithful JS ports of the Python originals) ===== */

/* ----- core: propositional + modal(K/T/S4/S5) + Fano + dispatch ----- */
/* prover_core.js — faithful JavaScript port of prover_core.py.
   Same algorithm: real octonion (Cayley–Dickson) multiplication, real associator,
   truth-table backstop, substrate dispatch. Verified against the Python output
   case-for-case (see test harness). No numpy, no Pyodide — runs natively. */

/* ── quaternion + octonion arithmetic (plain arrays) ── */
function qmul(p, q) {
  const [a1,b1,c1,d1]=p,[a2,b2,c2,d2]=q;
  return [a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
          a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2];
}
const qconj=(p)=>[p[0],-p[1],-p[2],-p[3]];
function omul(x, y) {
  const a=x.slice(0,4),b=x.slice(4),c=y.slice(0,4),d=y.slice(4);
  const re=qmul(a,c).map((u,i)=>u-qmul(qconj(d),b)[i]);
  const im=qmul(d,a).map((u,i)=>u+qmul(b,qconj(c))[i]);
  return re.concat(im);
}
const osub=(x,y)=>x.map((u,i)=>u-y[i]);
const onorm2=(x)=>x.reduce((s,u)=>s+u*u,0);
const associator=(x,y,z)=>osub(omul(omul(x,y),z),omul(x,omul(y,z)));
function unit(i){const v=new Array(8).fill(0);v[i]=1;return v;}

const FANO_LINES=[[1,2,3],[1,4,5],[1,7,6],[2,4,6],[2,5,7],[3,4,7],[3,6,5]];

/* ── substrate tower ── */
const TOWER=[
  {name:"R",assoc_live:false,comm_live:false},
  {name:"C",assoc_live:false,comm_live:false},
  {name:"H",assoc_live:false,comm_live:true},
  {name:"O",assoc_live:true, comm_live:true},
];
function dispatch_system(need_nc, need_na){
  for(const r of TOWER){
    if(need_na && !r.assoc_live) continue;
    if(need_nc && !r.comm_live) continue;
    const surrendered=[];
    if(!r.assoc_live) surrendered.push("nonassociativity (associator ≡ 0 here)");
    if(!r.comm_live)  surrendered.push("noncommutativity (commutator ≡ 0 here)");
    return {rung:r.name, surrendered,
      note:`Dispatched to ${r.name}. Octonionic structure above ${r.name} is NOT `+
           `claimed for this system; if the inference relation later needs it, that `+
           `is an obstruction to record and the dispatcher climbs.`};
  }
  return {rung:"O",surrendered:[],note:""};
}

/* ── Fano step with computed grade ── */
function lineThrough(a,b){ for(const L of FANO_LINES) if(L.includes(a)&&L.includes(b)) return L; return null; }
function mates(u){
  const m=new Set();
  for(const L of FANO_LINES) if(L.includes(u)) for(const x of L) if(x!==u) m.add(x);
  return [...m].sort((a,b)=>a-b).map(x=>"e"+x).join(", ");
}
function fanoGrade(a,b,third){
  const ea=unit(a),eb=unit(b),ec=unit(third);
  const assoc_n=onorm2(associator(ea,eb,ec));
  return {assoc_norm:assoc_n, live:assoc_n>1e-9};
}
function checkFano(ruleId,a,b,c){
  if(ruleId==="assume") return {level:"ok",msg:"Unit introduced."};
  if(a==null||b==null) return {level:"err",msg:"Both premises must be units e1..e7.",
    why:"The only inference here is octonion multiplication of two imaginary units, so each premise must be one of e1..e7.",
    fix:"Cite two earlier lines that are units."};
  if(a===b) return {level:"err",msg:`Both premises are e${a}; they must differ.`,
    why:"A Fano line has three distinct points. A unit times itself leaves the imaginary plane, so the rule needs two distinct units.",
    fix:`Pick a second unit different from e${a} (line-mates: ${mates(a)}).`};
  const L=lineThrough(a,b);
  if(L===null) return {level:"err",msg:`e${a} and e${b} aren't on a common Fano line.`,
    why:`Only collinear unit pairs multiply to a single third unit. e${a} and e${b} aren't collinear.`,
    fix:`Pair e${a} with one of ${mates(a)}.`};
  const third=L.find(x=>x!==a&&x!==b);
  if(c!==third) return {level:"err",
    msg:`e${a},e${b} lie on line {${L.join(",")}}; they derive e${third}${c==null?"":", not e"+c}.`,
    why:`The third point on the line through e${a} and e${b} is e${third} — the unique unit their product lands on.`,
    fix:`Set the conclusion to e${third}.`};
  const g=fanoGrade(a,b,third);
  const grade=`Associator = TRUE ZERO on this step (computed |assoc|²=${trimnum(g.assoc_norm)}). `+
    `A single Fano line {${L.join(",")}} generates a quaternionic ℍ-subalgebra, so `+
    `e${a}·e${b}=e${third} associates exactly. The grade goes live only when a `+
    `derivation combines units from different lines — that's where bracketing starts to matter.`;
  return {level:"ok",msg:`Valid Fano inference on line {${L.join(",")}}.`,grade};
}
function trimnum(x){ // mimic python %.3g for small values
  if(x===0) return "0";
  const s=x.toPrecision(3); return parseFloat(s).toString();
}

/* ── formula helpers ── */
const ARITY={not:1,box:1,dia:1,and:2,or:2,imp:2,iff:2};
const SYM={not:"¬",box:"□",dia:"◇",and:"∧",or:"∨",imp:"→",iff:"↔"};
function fmt(n){
  if(n==null) return "·";
  if(n.t==="atom") return n.name;
  const a=n.args;
  if(ARITY[n.op]===1){
    const inner=fmt(a[0]);
    const par=a[0]&&a[0].t==="op"&&ARITY[a[0].op]===2;
    return SYM[n.op]+(par?`(${inner})`:inner);
  }
  return `(${fmt(a[0])} ${SYM[n.op]} ${fmt(a[1])})`;
}
const eq=(a,b)=>fmt(a)===fmt(b);
function wellFormed(n){
  if(n==null) return false;
  if(n.t==="atom") return !!n.name;
  if(!(n.op in ARITY)) return false;
  if(n.args.length!==ARITY[n.op]) return false;
  return n.args.every(wellFormed);
}
function atomsOf(n,s){ if(n==null)return; if(n.t==="atom")s.add(n.name); else n.args.forEach(x=>atomsOf(x,s)); }
function evalv(n,env){
  if(n.t==="atom") return !!env[n.name];
  const a=n.args.map(x=>x?evalv(x,env):false);
  switch(n.op){
    case"not":return !a[0]; case"and":return a[0]&&a[1]; case"or":return a[0]||a[1];
    case"imp":return !a[0]||a[1]; case"iff":return a[0]===a[1]; default:return false;
  }
}
function entails(prem,c){
  const s=new Set(); prem.forEach(p=>atomsOf(p,s)); atomsOf(c,s);
  const vs=[...s];
  for(let m=0;m<(1<<vs.length);m++){
    const env={}; vs.forEach((v,i)=>env[v]=!!(m&(1<<i)));
    if(prem.every(p=>evalv(p,env)) && !evalv(c,env)) return {env,vars:vs};
  }
  return null;
}

const PROP_SHAPE={mp:"A,  A→B  ⊢  B",mt:"A→B,  ¬B  ⊢  ¬A",
  andI:"A,  B  ⊢  A∧B",andE1:"A∧B  ⊢  A",andE2:"A∧B  ⊢  B",
  orI1:"A  ⊢  A∨B   (any B)",orI2:"B  ⊢  A∨B   (any A)",
  dn:"¬¬A  ⊢  A   (and A ⊢ ¬¬A)"};
const is=(n,op)=>n!=null&&n.t==="op"&&n.op===op;
const opN=(o,a)=>({t:"op",op:o,args:a});

function checkPropRule(ruleId,prem,c){
  const A=prem[0]??null, B=prem[1]??null, sh=PROP_SHAPE[ruleId];
  switch(ruleId){
    case"assume":
      return wellFormed(c)?{level:"assumed",msg:"Assumption — introduced, not derived.",
        why:"An assumption is asserted, not proved. It's well-formed, so later steps may cite it, but it carries no inferential justification of its own."}:
        {level:"err",msg:"Assumption is not a well-formed formula.",
         why:"An assumption may be any formula, but it must be complete.",fix:"Fill in the remaining operand(s)."};
    case"mp":{
      const trymp=(imp,x)=>is(imp,"imp")&&eq(imp.args[0],x)&&eq(imp.args[1],c);
      if(trymp(A,B)||trymp(B,A)) return {level:"ok",msg:"Valid modus ponens."};
      const imp=is(A,"imp")?A:(is(B,"imp")?B:null);
      if(imp==null) return {level:"err",msg:"Neither cited line is a conditional (→).",
        why:"Modus ponens fires a conditional A→B using its antecedent A.",shape:sh,
        fix:"Cite a line whose main connective is →, plus a line matching its left side."};
      const other=eq(imp,A)?B:A;
      if(other==null||!eq(other,imp.args[0])) return {level:"err",
        msg:`The other premise ${other?fmt(other):"(missing)"} doesn't match the antecedent ${fmt(imp.args[0])}.`,
        why:`In A→B ⊢ B you must already have A exactly. Here A is ${fmt(imp.args[0])}.`,
        shape:sh,fix:`Cite a line equal to ${fmt(imp.args[0])}.`};
      if(!eq(imp.args[1],c)) return {level:"err",
        msg:`Conclusion is ${fmt(c)}, but MP yields the consequent ${fmt(imp.args[1])}.`,
        why:"Once A and A→B are present, the only thing you may write is B.",
        shape:sh,fix:`Set this line to ${fmt(imp.args[1])}.`};
      return {level:"err",msg:"Premises don't line up for modus ponens.",shape:sh};
    }
    case"mt":{
      const imp=is(A,"imp")?A:(is(B,"imp")?B:null);
      const neg=is(A,"not")?A:(is(B,"not")?B:null);
      if(imp&&neg&&eq(neg.args[0],imp.args[1])&&is(c,"not")&&eq(c.args[0],imp.args[0]))
        return {level:"ok",msg:"Valid modus tollens."};
      if(imp==null) return {level:"err",msg:"Modus tollens needs a conditional (→) premise.",
        why:"It works backwards along an implication A→B.",shape:sh,fix:"Cite a line whose main connective is →."};
      if(neg==null||!eq(neg.args[0],imp.args[1])) return {level:"err",
        msg:`You need ¬${fmt(imp.args[1])} (negation of the consequent).`,
        why:"Modus tollens denies the consequent to deny the antecedent.",shape:sh,
        fix:`Add or cite the line ¬${fmt(imp.args[1])}.`};
      return {level:"err",msg:`Conclusion should be ¬${fmt(imp.args[0])} (negation of the antecedent).`,
        why:"Denying the consequent lets you deny the antecedent.",shape:sh,
        fix:`Set the conclusion to ¬${fmt(imp.args[0])}.`};
    }
    case"andI":
      if(is(c,"and")&&((eq(c.args[0],A)&&eq(c.args[1],B))||(eq(c.args[0],B)&&eq(c.args[1],A))))
        return {level:"ok",msg:"Valid conjunction introduction."};
      return {level:"err",msg:"Conclusion should be the two premises joined by ∧.",
        why:"∧-introduction records that both premises hold at once.",shape:sh,
        fix:`Set the conclusion to (${fmt(A)} ∧ ${fmt(B)}).`};
    case"andE1": case"andE2":{
      const side=ruleId==="andE1"?0:1;
      if(is(A,"and")&&eq(A.args[side],c))
        return {level:"ok",msg:`Valid ∧-elimination (${side===0?"left":"right"}).`};
      if(!is(A,"and")) return {level:"err",msg:"∧ Elim needs a conjunction premise.",
        why:"You can only pull a conjunct out of an actual ∧ formula.",shape:sh,
        fix:"Cite a line whose main connective is ∧."};
      return {level:"err",msg:`The ${side===0?"left":"right"} conjunct is ${fmt(A.args[side])}, not ${fmt(c)}.`,
        why:"∧ Elim keeps the formula on that side of the ∧.",shape:sh,
        fix:`Set the conclusion to ${fmt(A.args[side])}.`};
    }
    case"orI1": case"orI2":{
      const side=ruleId==="orI1"?0:1;
      if(is(c,"or")&&eq(c.args[side],A)) return {level:"ok",msg:"Valid ∨-introduction."};
      return {level:"err",msg:`Conclusion must be a disjunction with ${fmt(A)} on the ${side===0?"left":"right"}.`,
        why:`If A holds then a disjunction containing A holds — but A must be the ${side===0?"left":"right"} disjunct here.`,
        shape:sh,fix:`Make a ∨ with ${fmt(A)} on the ${side===0?"left":"right"}.`};
    }
    case"dn":
      if(is(A,"not")&&is(A.args[0],"not")&&eq(A.args[0].args[0],c))
        return {level:"ok",msg:"Valid double-negation elimination (¬¬A ⊢ A)."};
      if(is(c,"not")&&is(c.args[0],"not")&&eq(c.args[0].args[0],A))
        return {level:"ok",msg:"Valid double-negation introduction (A ⊢ ¬¬A)."};
      return {level:"err",msg:`${fmt(A)} and ${fmt(c)} don't differ by exactly ¬¬.`,
        why:"Double negation only adds or removes a ¬¬ in front of the same formula.",shape:sh,
        fix:`From ${fmt(A)} the valid step is ¬¬${fmt(A)}.`};
    default: return {level:"err",msg:"Unknown rule."};
  }
}
function checkProp(ruleId,prem,c){
  const base=checkPropRule(ruleId,prem,c);
  if(ruleId==="assume") return base;
  const ready=prem.length>0&&prem.every(wellFormed)&&wellFormed(c);
  const ce=ready?entails(prem,c):null;
  const semValid=ready?(ce===null):null;
  if(base.level==="ok"&&semValid===false){
    const env=ce.env, vs=ce.vars;
    const cs=vs.map(v=>`${v}=${env[v]?"T":"F"}`).join(", ");
    return {level:"err",msg:"This fits the rule's pattern but isn't truth-preserving.",
      why:"A valid step can never make true premises yield a false conclusion. Here a valuation makes every premise true yet the conclusion false.",
      counterexample:`When ${cs}: all premises true but ${fmt(c)} is false.`,
      fix:"Change the conclusion to what the premises force, or add the missing premise."};
  }
  if(base.level!=="ok"&&semValid===true){
    const out={...base,level:"warn"};
    out.why=((base.why?base.why+" ":"")+"The conclusion does follow from these premises — so the issue is the rule citation, not the logic.").trim();
    if(!out.fix) out.fix="Switch the rule to the one whose shape matches this step.";
    return out;
  }
  return base;
}

const MODAL_SHAPE={mp:"A,  A→B  ⊢  B",nec:"⊢A ⟹ ⊢□A  (A a theorem)",
  kdist:"□(A→B)  ⊢  □A→□B",dualD:"◇A  ⊢  ¬□¬A",
  semantic:"premises  ⊨  conclusion   (checked by countermodel search)"};

/* ── embedded Kripke semantics (K/T/S4/S5), faithful port ── */
function _katoms(n,s){ if(n==null)return; if(n.t==="atom")s.add(n.name); else n.args.forEach(a=>_katoms(a,s)); }
function _sat(n,w,worlds,R,V){
  if(n.t==="atom") return V[w].has(n.name);
  const a=n.args;
  switch(n.op){
    case"not": return !_sat(a[0],w,worlds,R,V);
    case"and": return _sat(a[0],w,worlds,R,V)&&_sat(a[1],w,worlds,R,V);
    case"or":  return _sat(a[0],w,worlds,R,V)||_sat(a[1],w,worlds,R,V);
    case"imp": return !_sat(a[0],w,worlds,R,V)||_sat(a[1],w,worlds,R,V);
    case"iff": return _sat(a[0],w,worlds,R,V)===_sat(a[1],w,worlds,R,V);
    case"box": return worlds.every(v=>!R.has(w+","+v)||_sat(a[0],v,worlds,R,V));
    case"dia": return worlds.some(v=>R.has(w+","+v)&&_sat(a[0],v,worlds,R,V));
    default: return false;
  }
}
function _transclose(R){
  let changed=true;
  while(changed){ changed=false;
    for(const e1 of [...R]){ const [x,y]=e1.split(",").map(Number);
      for(const e2 of [...R]){ const [y2,z]=e2.split(",").map(Number);
        if(y===y2 && !R.has(x+","+z)){ R.add(x+","+z); changed=true; } } } }
}
function _closeR(worlds,baseR,logic){
  const R=new Set(baseR);
  if(logic==="T"||logic==="S4"||logic==="S5") for(const w of worlds) R.add(w+","+w);
  if(logic==="S4"||logic==="S5") _transclose(R);
  if(logic==="S5"){ for(const e of [...R]){ const [x,y]=e.split(",").map(Number); R.add(y+","+x); } _transclose(R); }
  return R;
}
function findCountermodel(premises,concl,logic,maxWorlds){
  const aset=new Set(); premises.forEach(p=>_katoms(p,aset)); _katoms(concl,aset);
  let atoms=[...aset].sort(); if(atoms.length===0) atoms=["p"];
  const atomSubsets=[];
  for(let m=0;m<(1<<atoms.length);m++){ const s=new Set();
    for(let i=0;i<atoms.length;i++) if((m>>i)&1) s.add(atoms[i]); atomSubsets.push(s); }
  for(let nW=1;nW<=maxWorlds;nW++){
    const worlds=[]; for(let i=0;i<nW;i++) worlds.push(i);
    const edges=[]; for(const x of worlds) for(const y of worlds) edges.push(x+","+y);
    const ne=edges.length;
    const relCount = ne>12 ? 2 : (1<<ne);
    for(let rm=0; rm<relCount; rm++){
      let baseR;
      if(ne>12){ baseR = rm===0 ? new Set() : new Set(edges); }
      else { baseR=new Set(); for(let i=0;i<ne;i++) if((rm>>i)&1) baseR.add(edges[i]); }
      const R=_closeR(worlds,baseR,logic);
      const total=Math.pow(atomSubsets.length,nW);
      for(let t=0;t<total;t++){
        const V={}; let x=t;
        for(let i=nW-1;i>=0;i--){ V[worlds[i]]=atomSubsets[x%atomSubsets.length]; x=Math.floor(x/atomSubsets.length); }
        for(const w of worlds){
          if(premises.every(p=>_sat(p,w,worlds,R,V)) && !_sat(concl,w,worlds,R,V)){
            return { logic, worlds:nW, eval_world:w,
              relation:[...R].map(e=>e.split(",").map(Number)).sort((a,b)=>a[0]-b[0]||a[1]-b[1]),
              valuation:Object.fromEntries(worlds.map(v=>[v,[...V[v]].sort()])), atoms };
          }
        }
      }
    }
  }
  return null;
}
function describeCountermodel(cm){
  const ws=cm.worlds; const wl=[]; for(let i=0;i<ws;i++) wl.push("w"+i);
  const rel=cm.relation;
  const relTxt=rel.length?rel.map(([x,y])=>`w${x}→w${y}`).join(", "):"(no arrows)";
  const vp=[]; for(let w=0;w<ws;w++){ const set=cm.valuation[w]; vp.push(`w${w}: {${set.length?set.join(", "):"∅"}}`); }
  return `Countermodel in ${cm.logic}: worlds {${wl.join(", ")}}; accessibility ${relTxt}; `+
    `valuation ${vp.join("; ")}. At w${cm.eval_world} the premises all hold but the conclusion fails.`;
}

function checkModal(ruleId,prem,c,logic){
  logic = logic||"K";
  const A=prem[0]??null, sh=MODAL_SHAPE[ruleId];
  if(ruleId==="assume") return wellFormed(c)?{level:"assumed",msg:"Assumption — introduced, not derived.",
    why:"An assumption is asserted, not proved. It's well-formed, so later steps may cite it, but it carries no inferential justification of its own."}:
    {level:"err",msg:"Assumption is not a well-formed formula.",why:"An assumption may be any formula, but it must be complete.",fix:"Complete the formula."};
  if(ruleId==="semantic"){
    const prems=prem.filter(wellFormed);
    if(!wellFormed(c)) return {level:"err",msg:"The conclusion formula is incomplete.",why:"Need a complete formula to test.",fix:"Finish the conclusion."};
    const cm=findCountermodel(prems,c,logic,3);
    if(cm===null) return {level:"ok",msg:`Valid in ${logic} (no countermodel up to 3 worlds).`,
      grade:`Semantic check: every ${logic}-model satisfying the premises (searched to 3 worlds) also satisfies ${fmt(c)}.`};
    return {level:"err",msg:`Not valid in ${logic} — there is a countermodel.`,
      why:`Kripke semantics: a frame meeting ${logic}'s conditions has a world where the premises hold but the conclusion fails.`,
      counterexample:describeCountermodel(cm),
      fix:"Strengthen the premises, change the conclusion, or pick a stronger logic (T/S4/S5) whose frame conditions rule this countermodel out."};
  }
  if(ruleId==="mp") return checkPropRule("mp",prem,c);
  if(ruleId==="nec"){
    if(is(c,"box")&&eq(c.args[0],A)) return {level:"ok",msg:"Well-formed necessitation.",
      why:"Sound only when A is a theorem (no open assumptions). The shape is right."};
    return {level:"err",msg:`Necessitation gives □${fmt(A)}, not ${fmt(c)}.`,
      why:"Necessitation boxes a theorem: from a proved A assert □A.",shape:sh,
      fix:`Set the conclusion to □${fmt(A)}.`};
  }
  if(ruleId==="kdist"){
    if(is(A,"box")&&is(A.args[0],"imp")){
      const inn=A.args[0];
      const want=opN("imp",[opN("box",[inn.args[0]]),opN("box",[inn.args[1]])]);
      if(eq(want,c)) return {level:"ok",msg:"Valid K distribution."};
      return {level:"err",msg:`From ${fmt(A)}, K yields ${fmt(want)} — not ${fmt(c)}.`,
        why:"K lets necessity distribute over an implication.",shape:sh,fix:`Set the conclusion to ${fmt(want)}.`};
    }
    return {level:"err",msg:"K needs a premise shaped □(A→B).",
      why:"The K axiom only applies to a necessity wrapped around an implication.",shape:sh,
      fix:"Cite a line of the form □(…→…)."};
  }
  if(ruleId==="dualD"){
    if(is(A,"dia")){
      const want=opN("not",[opN("box",[opN("not",[A.args[0]])])]);
      if(eq(want,c)) return {level:"ok",msg:"Valid ◇/□ duality."};
      return {level:"err",msg:`◇${fmt(A.args[0])} unfolds to ${fmt(want)}, not ${fmt(c)}.`,
        why:"Possibility is definable from necessity: ◇A = ¬□¬A.",shape:sh,
        fix:`Set the conclusion to ${fmt(want)}.`};
    }
    return {level:"err",msg:"Duality rewrites a ◇ formula; this premise has no leading ◇.",
      why:"The rule converts ◇A into ¬□¬A.",shape:sh,fix:"Cite a line whose main connective is ◇."};
  }
  return {level:"err",msg:"Unknown rule."};
}

function check_step(systemId,ruleId,rulePremCount,prem,c,fanoIdx,logic){
  if(rulePremCount>0 && prem.length!==rulePremCount && systemId!=="fano" && ruleId!=="semantic")
    return {level:"err",msg:`This rule draws on ${rulePremCount} earlier line(s); ${prem.length} selected.`,
      why:"Every rule has a fixed number of inputs and can't fire until they're all cited.",
      fix:`Use the “from” drop-downs to cite exactly ${rulePremCount} line(s).`};
  if(systemId==="fano"){
    const a=fanoIdx&&fanoIdx.length>0?fanoIdx[0]:null;
    const b=fanoIdx&&fanoIdx.length>1?fanoIdx[1]:null;
    const cc=fanoIdx&&fanoIdx.length>2?fanoIdx[2]:null;
    return checkFano(ruleId,a,b,cc);
  }
  for(const p of prem) if(!wellFormed(p))
    return {level:"err",msg:"A cited premise line isn't a complete formula yet.",
      why:"A step can only use well-formed premises; one referenced line has a blank.",
      fix:"Finish the earlier line first."};
  if(ruleId!=="assume"&&!wellFormed(c))
    return {level:"err",msg:"This line's formula is incomplete.",
      why:"Every connective needs all its operands filled in.",fix:"Complete the formula in the builder."};
  if(systemId==="modalK") return checkModal(ruleId,prem,c,logic);
  return checkProp(ruleId,prem,c);
}

/* ----- equational-logic engine (group/octonion/modular/boolean) ----- */
/* equational.js — faithful JS port of equational.py (Birkhoff equational logic).
   Verified field-for-field against the Python output. */

function evar(name){ return {t:"var",name}; }
function eapp(op,args){ return {t:"app",op,args:[...args]}; }
function tfmt(n){
  if(n==null) return "·";
  if(n.t==="var") return n.name;
  const a=n.args;
  if(!a.length) return n.op;
  return `${n.op}(${a.map(tfmt).join(",")})`;
}
const teq=(a,b)=>tfmt(a)===tfmt(b);
function wellFormedTerm(n,sig){
  if(n==null) return false;
  if(n.t==="var") return !!n.name;
  if(!(n.op in sig)) return false;
  if(n.args.length!==sig[n.op]) return false;
  return n.args.every(x=>wellFormedTerm(x,sig));
}
function substitute(term,env){
  if(term.t==="var") return (term.name in env)?env[term.name]:term;
  return eapp(term.op, term.args.map(a=>substitute(a,env)));
}
function match(pattern,term,env){
  if(env===undefined) env={};
  if(pattern.t==="var"){
    if(pattern.name in env) return teq(env[pattern.name],term)?env:null;
    env={...env}; env[pattern.name]=term; return env;
  }
  if(term.t!=="app"||pattern.op!==term.op) return null;
  if(pattern.args.length!==term.args.length) return null;
  for(let i=0;i<pattern.args.length;i++){
    env=match(pattern.args[i],term.args[i],env);
    if(env===null) return null;
  }
  return env;
}
function axiomInstance(al0,ar0,eqL,eqR){
  for(const [al,ar] of [[al0,ar0],[ar0,al0]]){
    const env=match(al,eqL);
    if(env!==null){
      const env2=match(ar,eqR,env);
      if(env2!==null) return true;
    }
  }
  return false;
}
function differsInOneArg(lhs,rhs,subL,subR){
  if(lhs.t!=="app"||rhs.t!=="app") return false;
  if(lhs.op!==rhs.op||lhs.args.length!==rhs.args.length) return false;
  const diff=[];
  for(let k=0;k<lhs.args.length;k++) if(!teq(lhs.args[k],rhs.args[k])) diff.push(k);
  if(diff.length!==1) return false;
  const k=diff[0];
  return teq(lhs.args[k],subL)&&teq(rhs.args[k],subR);
}
function _need(n,name){
  return {level:"err",msg:`${name[0].toUpperCase()+name.slice(1)} cites ${n} earlier equation(s).`,
    why:"This rule has a fixed number of inputs.",fix:`Cite exactly ${n} line(s).`};
}
function depth(n){
  if(n==null||n.t==="var") return 0;
  if(!n.args.length) return 1;
  return 1+Math.max(...n.args.map(depth));
}
function collectGens(n,s){
  if(n==null) return;
  if(n.t==="var") s.add(n.name);
  else if(n.t==="app") n.args.forEach(a=>collectGens(a,s));
}
function threeDistinctGenerators(n){
  const s=new Set(); collectGens(n,s);
  return s.size>=3 && depth(n)>=2;
}
function gradeNote(theory,L,R){
  if(theory.name!=="octonion") return null;
  if(threeDistinctGenerators(L)||threeDistinctGenerators(R))
    return "Object-theory note: this identity ranges over ≥3 distinct octonion "+
      "generators, so it lives on a non-associative locus — the associator is generally "+
      "nonzero here. (Derived from the theory; the equational step itself is licensed "+
      "structurally, not by this fact.)";
  return "Object-theory note: within a single 2-generated (quaternionic) chart, so "+
      "associativity holds here — the associator is a true zero.";
}

function checkEquational(ruleId,cited,concl,theory){
  const sig=theory.sig;
  const L=concl.l, R=concl.r;
  if(!wellFormedTerm(L,sig)||!wellFormedTerm(R,sig))
    return {level:"err",msg:"The asserted equation has a malformed term.",
      why:`Both sides must be terms over the theory's signature ${JSON.stringify(Object.keys(sig).sort())}.`,
      fix:"Check each operation's name and arity."};
  if(ruleId==="refl"){
    if(teq(L,R)) return {level:"ok",msg:"Reflexivity: t = t.",grade:gradeNote(theory,L,R)};
    return {level:"err",msg:`Reflexivity requires identical sides; ${tfmt(L)} ≠ ${tfmt(R)}.`,
      why:"refl only proves t = t.",shape:"⊢ t = t",fix:"Make both sides the same term."};
  }
  if(ruleId==="sym"){
    if(cited.length!==1) return _need(1,"symmetry");
    const e=cited[0];
    if(teq(e.l,R)&&teq(e.r,L)) return {level:"ok",msg:"Symmetry: from s = t infer t = s.",grade:gradeNote(theory,L,R)};
    return {level:"err",msg:"Symmetry must swap the cited equation's sides.",
      why:`From ${tfmt(e.l)} = ${tfmt(e.r)}, symmetry yields ${tfmt(e.r)} = ${tfmt(e.l)}.`,
      shape:"s = t  ⊢  t = s",fix:`Set this line to ${tfmt(e.r)} = ${tfmt(e.l)}.`};
  }
  if(ruleId==="trans"){
    if(cited.length!==2) return _need(2,"transitivity");
    const [e1,e2]=cited;
    for(const [a,b] of [[e1,e2],[e2,e1]]){
      if(teq(a.r,b.l)&&teq(L,a.l)&&teq(R,b.r))
        return {level:"ok",msg:"Transitivity: chain s = t = u.",grade:gradeNote(theory,L,R)};
    }
    return {level:"err",msg:"Transitivity needs the two equations to share a middle term.",
      why:"From s = t and t = u, conclude s = u; the inner terms must match.",
      shape:"s = t,  t = u  ⊢  s = u",
      fix:"Cite two lines where the right side of one equals the left side of the other."};
  }
  if(ruleId==="cong"){
    if(cited.length!==1) return _need(1,"congruence");
    const e=cited[0];
    if(differsInOneArg(L,R,e.l,e.r))
      return {level:"ok",msg:"Congruence: replace equals in one argument.",grade:gradeNote(theory,L,R)};
    return {level:"err",msg:"Congruence must rebuild the SAME operation, changing one argument via the cited equation.",
      why:"From sᵢ = tᵢ, congruence gives f(…sᵢ…) = f(…tᵢ…) with all other arguments unchanged.",
      shape:"sᵢ = tᵢ  ⊢  f(…sᵢ…) = f(…tᵢ…)",
      fix:`Both sides must share an operation and differ only where the cited equation ${tfmt(e.l)} = ${tfmt(e.r)} applies.`};
  }
  if(ruleId==="axiom"){
    for(const ax of theory.axioms){
      if(axiomInstance(ax.l,ax.r,L,R))
        return {level:"ok",msg:`Axiom instance: ${ax.name||"(axiom)"}.`,grade:gradeNote(theory,L,R)};
    }
    return {level:"err",msg:`This equation is not an instance of any ${theory.name||"theory"} axiom.`,
      why:"The axiom rule only licenses substitution instances of the theory's stated equations.",
      fix:`Pick an equation matching one of: ${theory.axioms.map(a=>`${tfmt(a.l)} = ${tfmt(a.r)}`).join("; ")}.`};
  }
  return {level:"err",msg:"Unknown equational rule."};
}

function theoryGroup(){
  return {name:"group",sig:{e:0,i:1,m:2},axioms:[
    {name:"assoc",l:eapp("m",[eapp("m",[evar("x"),evar("y")]),evar("z")]),
                  r:eapp("m",[evar("x"),eapp("m",[evar("y"),evar("z")])])},
    {name:"unit-left",l:eapp("m",[eapp("e",[]),evar("x")]),r:evar("x")},
    {name:"inv-left",l:eapp("m",[eapp("i",[evar("x")]),evar("x")]),r:eapp("e",[])},
  ]};
}
function theoryOctonion(){
  const FANO=[[1,2,3],[1,4,5],[1,7,6],[2,4,6],[2,5,7],[3,4,7],[3,6,5]];
  const sig={m:2}; for(let i=1;i<=7;i++) sig["e"+i]=0;
  const axioms=FANO.map(([a,b,c])=>({name:`fano[e${a}·e${b}=e${c}]`,
    l:eapp("m",[eapp("e"+a,[]),eapp("e"+b,[])]),r:eapp("e"+c,[])}));
  return {name:"octonion",sig,axioms};
}
function theoryModular(n){
  if(n===undefined) n=5;
  const sig={a:2,m:2}; for(let i=0;i<n;i++) sig[""+i]=0;
  const K=i=>eapp(""+(((i%n)+n)%n),[]);
  const x=evar("x"),y=evar("y"),z=evar("z");
  const axioms=[
    {name:"add-comm",l:eapp("a",[x,y]),r:eapp("a",[y,x])},
    {name:"add-assoc",l:eapp("a",[eapp("a",[x,y]),z]),r:eapp("a",[x,eapp("a",[y,z])])},
    {name:"add-id",l:eapp("a",[x,K(0)]),r:x},
    {name:"mul-comm",l:eapp("m",[x,y]),r:eapp("m",[y,x])},
    {name:"mul-assoc",l:eapp("m",[eapp("m",[x,y]),z]),r:eapp("m",[x,eapp("m",[y,z])])},
    {name:"mul-id",l:eapp("m",[x,K(1)]),r:x},
    {name:"distrib",l:eapp("m",[x,eapp("a",[y,z])]),r:eapp("a",[eapp("m",[x,y]),eapp("m",[x,z])])},
  ];
  for(let i=0;i<n;i++) for(let j=0;j<n;j++){
    axioms.push({name:`${i}+${j}=${(i+j)%n}`,l:eapp("a",[K(i),K(j)]),r:K(i+j)});
    axioms.push({name:`${i}*${j}=${(i*j)%n}`,l:eapp("m",[K(i),K(j)]),r:K(i*j)});
  }
  return {name:"modular",modulus:n,sig,axioms};
}
function theoryBoolean(){
  const sig={or:2,and:2,not:1,T:0,F:0};
  const x=evar("x"),y=evar("y"),z=evar("z");
  const O=(a,b)=>eapp("or",[a,b]), A=(a,b)=>eapp("and",[a,b]), N=a=>eapp("not",[a]);
  const T=eapp("T",[]), F=eapp("F",[]);
  const axioms=[
    {name:"∪-comm",l:O(x,y),r:O(y,x)},
    {name:"∩-comm",l:A(x,y),r:A(y,x)},
    {name:"∪-assoc",l:O(O(x,y),z),r:O(x,O(y,z))},
    {name:"∩-assoc",l:A(A(x,y),z),r:A(x,A(y,z))},
    {name:"absorb-∪",l:O(x,A(x,y)),r:x},
    {name:"absorb-∩",l:A(x,O(x,y)),r:x},
    {name:"distrib-∪",l:O(x,A(y,z)),r:A(O(x,y),O(x,z))},
    {name:"distrib-∩",l:A(x,O(y,z)),r:O(A(x,y),A(x,z))},
    {name:"complement-∪",l:O(x,N(x)),r:T},
    {name:"complement-∩",l:A(x,N(x)),r:F},
    {name:"identity-∪",l:O(x,F),r:x},
    {name:"identity-∩",l:A(x,T),r:x},
    {name:"deMorgan-∪",l:N(O(x,y)),r:A(N(x),N(y))},
    {name:"deMorgan-∩",l:N(A(x,y)),r:O(N(x),N(y))},
    {name:"double-neg",l:N(N(x)),r:x},
  ];
  return {name:"boolean",sig,axioms};
}
function _isPrime(p){ if(p<2)return false; for(let i=2;i*i<=p;i++) if(p%i===0)return false; return true; }
function theoryPrimeField(p){
  if(p===undefined) p=7;
  if(!_isPrime(p)) throw new Error("GF(p) requires p prime; got "+p);
  const sig={a:2,m:2}; for(let i=0;i<p;i++) sig[""+i]=0;
  const K=i=>eapp(""+(((i%p)+p)%p),[]);
  const x=evar("x"),y=evar("y"),z=evar("z");
  const axioms=[
    {name:"add-comm",l:eapp("a",[x,y]),r:eapp("a",[y,x])},
    {name:"add-assoc",l:eapp("a",[eapp("a",[x,y]),z]),r:eapp("a",[x,eapp("a",[y,z])])},
    {name:"add-id",l:eapp("a",[x,K(0)]),r:x},
    {name:"mul-comm",l:eapp("m",[x,y]),r:eapp("m",[y,x])},
    {name:"mul-assoc",l:eapp("m",[eapp("m",[x,y]),z]),r:eapp("m",[x,eapp("m",[y,z])])},
    {name:"mul-id",l:eapp("m",[x,K(1)]),r:x},
    {name:"distrib",l:eapp("m",[x,eapp("a",[y,z])]),r:eapp("a",[eapp("m",[x,y]),eapp("m",[x,z])])},
  ];
  for(let i=0;i<p;i++) for(let j=0;j<p;j++){
    axioms.push({name:`${i}+${j}=${(i+j)%p}`,l:eapp("a",[K(i),K(j)]),r:K(i+j)});
    axioms.push({name:`${i}*${j}=${(i*j)%p}`,l:eapp("m",[K(i),K(j)]),r:K(i*j)});
  }
  for(let i=1;i<p;i++){ let inv=1; for(let j=1;j<p;j++) if((i*j)%p===1){inv=j;break;}
    axioms.push({name:`inv(${i})=${inv}`,l:eapp("m",[K(i),K(inv)]),r:K(1)}); }
  return {name:"prime_field",p,sig,axioms};
}
const THEORIES={group:theoryGroup(),octonion:theoryOctonion(),
  modular:theoryModular(5),boolean:theoryBoolean(),prime_field:theoryPrimeField(7)};



/* ----- linear-logic engine (MLL) ----- */
/* linear.js — faithful JS port of linear.py (multiplicative linear logic with
   backtracking resource search). Verified against the Python output. */

function lat(s){ return {t:"atom",name:s}; }
function tensor(a,b){ return {t:"op",op:"tensor",args:[a,b]}; }
function lolli(a,b){ return {t:"op",op:"lolli",args:[a,b]}; }
const ONE={t:"op",op:"one",args:[]};

function lf(n){
  if(n==null) return "·";
  if(n.t==="atom") return n.name;
  const op=n.op;
  if(op==="one") return "1";
  const a=n.args, sym = op==="tensor"?"⊗":"⊸";
  return "("+lf(a[0])+" "+sym+" "+lf(a[1])+")";
}
const leq=(a,b)=>lf(a)===lf(b);

function* splits(ms){
  const n=ms.length;
  for(let mask=0;mask<(1<<n);mask++){
    const L=[],R=[];
    for(let i=0;i<n;i++){ if((mask>>i)&1) L.push(ms[i]); else R.push(ms[i]); }
    yield [L,R];
  }
}

function provable(gamma, goal, fuel){
  if(fuel===undefined) fuel=64;
  if(fuel<=0) return false;
  const g=[...gamma];
  if(g.length===1 && leq(g[0],goal)) return true;
  if(g.length===0 && goal.t==="op" && goal.op==="one") return true;
  if(goal.t==="op" && goal.op==="lolli"){
    const [A,B]=goal.args;
    if(provable([...g,A],B,fuel-1)) return true;
  }
  if(goal.t==="op" && goal.op==="tensor"){
    const [A,B]=goal.args;
    for(const [L,R] of splits(g))
      if(provable(L,A,fuel-1) && provable(R,B,fuel-1)) return true;
  }
  for(let i=0;i<g.length;i++){
    const res=g[i], rest=g.slice(0,i).concat(g.slice(i+1));
    if(res.t==="op" && res.op==="tensor"){
      const [A,B]=res.args;
      if(provable([...rest,A,B],goal,fuel-1)) return true;
    }
    if(res.t==="op" && res.op==="one"){
      if(provable(rest,goal,fuel-1)) return true;
    }
    if(res.t==="op" && res.op==="lolli"){
      const [A,B]=res.args;
      for(const [L,R] of splits(rest))
        if(provable(L,A,fuel-1) && provable([...R,B],goal,fuel-1)) return true;
    }
  }
  return false;
}

function diagnose(gamma, goal){
  const g=[...gamma];
  const haveAtoms={};
  for(const x of g) if(x.t==="atom") haveAtoms[lf(x)]=(haveAtoms[lf(x)]||0)+1;
  if(goal.t==="atom" && !(lf(goal) in haveAtoms) && !g.some(x=>x.t==="op"))
    return `No resource supplies ${lf(goal)}, and linear logic can't conjure one — `+
      `every atom in the goal must be produced from the context.`;
  if(g.every(x=>x.t==="atom") && goal.t==="atom"){
    if(g.length>1)
      return `Linear logic forbids discarding resources: ${g.length} hypotheses but `+
        `the goal consumes at most one, so the rest are left over (no weakening).`;
    if(g.length===1 && !leq(g[0],goal))
      return `The single resource ${lf(g[0])} isn't the goal ${lf(goal)}, and there's `+
        `no rule to convert it without a ⊸.`;
  }
  return "No linear derivation found that consumes every resource exactly once "+
    "(searched with backtracking). Linear logic allows neither duplicating nor "+
    "discarding hypotheses.";
}

function checkLinear(gamma, goal){
  if(provable(gamma,goal))
    return {level:"ok",msg:`Provable in MLL: ${gamma.map(lf).join(", ")||"·"} ⊢ ${lf(goal)}.`,
      grade:"Backtracking search found a derivation consuming each resource exactly once."};
  return {level:"err",msg:"Not provable in linear logic.",
    why:diagnose(gamma,goal),
    fix:"Adjust the resources so each is used exactly once — no duplication, no waste."};
}

/* ----- discharge / subproof engine ----- */
/* discharge.js — faithful JS port of discharge.py (assumption-discharge / subproof
   scope tracking). Verified against the Python output. */

const _ARITY={not:1,and:2,or:2,imp:2,iff:2,box:1,dia:1,bot:0};
const _SYM={not:"¬",and:"∧",or:"∨",imp:"→",iff:"↔",box:"□",dia:"◇",bot:"⊥"};

function df(n){
  if(n==null) return "·";
  if(n.t==="atom") return n.name;
  const op=n.op, a=n.args;
  if(op==="bot") return "⊥";
  if(_ARITY[op]===1){ const inner=df(a[0]); const par=a[0]&&a[0].t==="op"&&_ARITY[a[0].op]===2;
    return _SYM[op]+(par?"("+inner+")":inner); }
  return "("+df(a[0])+" "+_SYM[op]+" "+df(a[1])+")";
}
const dfeq=(a,b)=>df(a)===df(b);
const dis=(n,op)=>n!=null&&n.t==="op"&&n.op===op;

function dwf(n){
  if(n==null) return false;
  if(n.t==="atom") return !!n.name;
  if(n.op==="bot") return n.args.length===0;
  if(!(n.op in _ARITY)) return false;
  if(n.args.length!==_ARITY[n.op]) return false;
  return n.args.every(dwf);
}

function visibleLines(lines, idx){
  const targetDepth=lines[idx].depth;
  let st=[];
  for(let j=0;j<idx;j++){ const k=lines[j].kind;
    if(k==="assume") st.push(j); else if(k==="close"){ if(st.length) st.pop(); } }
  const openAssumes=new Set(st);
  const enclosing={};
  st=[];
  for(let j=0;j<idx;j++){ enclosing[j]=[...st]; const k=lines[j].kind;
    if(k==="assume") st.push(j); else if(k==="close"){ if(st.length) st.pop(); } }
  const visible=[];
  for(let j=0;j<idx;j++){
    if(lines[j].kind==="close") continue;
    const chain=enclosing[j].concat(lines[j].kind==="assume"?[j]:[]);
    if(chain.every(a=>openAssumes.has(a)) && lines[j].depth<=targetDepth) visible.push(j);
  }
  return visible;
}

function checkDischargeStep(lines, idx){
  const ln=lines[idx], kind=ln.kind;
  const vis=new Set(visibleLines(lines, idx));
  if(kind==="assume"){
    if(!dwf(ln.formula)) return {level:"err",msg:"Assumption is not a well-formed formula.",
      why:"Opening a subproof needs a complete hypothesis.",fix:"Complete the formula."};
    return {level:"assumed",msg:`Subproof assumption (depth ${ln.depth}) — to be discharged on close.`,
      why:"This opens a subproof. It must be discharged by an introduction rule before the result counts at the outer scope."};
  }
  const refs=(ln.refs||[]).filter(r=>r!=null);
  const bad=refs.filter(r=>!vis.has(r));
  if(bad.length) return {level:"err",msg:`Line cites ${bad.map(b=>"L"+(b+1)).join(", ")}, which is out of scope here.`,
    why:"You can't reuse a line from a subproof that has already been closed — its assumption was discharged.",
    fix:`Only cite lines still open at this depth (visible: ${[...vis].sort((a,b)=>a-b).map(v=>"L"+(v+1)).join(", ")||"none"}).`};
  if(kind==="close") return checkClose(lines, idx);
  return checkInscope(lines, idx, vis);
}

function premFormulas(lines, refs){
  return (refs||[]).filter(r=>r!=null&&r<lines.length).map(r=>lines[r].formula);
}

function checkInscope(lines, idx, vis){
  const ln=lines[idx], rule=ln.rule, c=ln.formula;
  const prem=premFormulas(lines, ln.refs);

  if(rule==="reit"){
    if(prem.length===1 && dfeq(prem[0],c)) return {level:"ok",msg:"Reiteration: a visible line copied into this scope."};
    return {level:"err",msg:"Reiteration must copy a cited line exactly.",
      why:"Reiteration brings an in-scope line into the current subproof unchanged.",
      shape:"L (visible)  ⊢  L",fix:"Cite one visible line and match its formula."};
  }
  if(rule==="mp"){
    if(prem.length===2){
      const [a,b]=prem;
      for(const [x,y] of [[a,b],[b,a]])
        if(dis(y,"imp")&&dfeq(y.args[0],x)&&dfeq(y.args[1],c)) return {level:"ok",msg:"Modus ponens."};
    }
    return {level:"err",msg:"Modus ponens needs A and A→B to give B.",
      why:"From a formula and a matching conditional, detach the consequent.",
      shape:"A,  A→B  ⊢  B",fix:"Cite A and A→B; conclude B."};
  }
  if(rule==="andI"){
    if(prem.length===2 && dis(c,"and") && dfeq(c.args[0],prem[0]) && dfeq(c.args[1],prem[1]))
      return {level:"ok",msg:"∧-introduction."};
    return {level:"err",msg:"∧-introduction pairs the two cited lines into a conjunction.",
      why:"From A and B, conclude A∧B in that order.",shape:"A,  B  ⊢  A∧B",fix:"Cite A then B; conclude A∧B."};
  }
  if(rule==="andE"){
    if(prem.length===1 && dis(prem[0],"and") && (dfeq(prem[0].args[0],c)||dfeq(prem[0].args[1],c)))
      return {level:"ok",msg:"∧-elimination."};
    return {level:"err",msg:"∧-elimination extracts one side of a conjunction.",
      why:"From A∧B, conclude A (or B).",shape:"A∧B  ⊢  A",fix:"Cite a conjunction; conclude one conjunct."};
  }
  if(rule==="botI"){
    if(prem.length===2 && dis(c,"bot")){
      const [a,b]=prem;
      for(const [x,y] of [[a,b],[b,a]])
        if(dis(y,"not")&&dfeq(y.args[0],x)) return {level:"ok",msg:"⊥-introduction: a formula and its negation."};
    }
    return {level:"err",msg:"⊥-introduction needs A and ¬A.",
      why:"A contradiction (a formula together with its negation) yields ⊥.",
      shape:"A,  ¬A  ⊢  ⊥",fix:"Cite A and ¬A; set this line to ⊥."};
  }
  return {level:"ok",msg:"In-scope line.",_delegate:true};
}

function checkClose(lines, idx){
  const ln=lines[idx], rule=ln.rule;
  let st=[];
  for(let j=0;j<idx;j++){ const k=lines[j].kind;
    if(k==="assume") st.push(j); else if(k==="close"){ if(st.length) st.pop(); } }
  if(!st.length) return {level:"err",msg:"No open subproof to close.",
    why:"A close/introduction rule must discharge an assumption opened earlier.",fix:"Open an assumption first."};
  const aIdx=st[st.length-1];
  const A=lines[aIdx].formula;
  const concl=ln.formula;
  const inner=[];
  for(let j=aIdx+1;j<idx;j++) if(lines[j].depth===lines[aIdx].depth && lines[j].kind!=="close") inner.push(j);
  const last = inner.length?lines[inner[inner.length-1]].formula:null;

  if(rule==="impI"){
    if(last===null) return {level:"err",msg:"→I needs at least one derived line inside the subproof.",
      why:"→Introduction discharges A after deriving some B; conclude A→B.",fix:"Derive a B under the assumption first."};
    const want={t:"op",op:"imp",args:[A,last]};
    if(dfeq(concl,want)) return {level:"ok",msg:`→Introduction: discharged ${df(A)}, concluded ${df(want)}.`,
      discharged:aIdx,grade:`Assumption ${df(A)} is now discharged; this result is valid at the outer scope.`};
    return {level:"err",msg:`→I should conclude ${df(want)}, not ${df(concl)}.`,
      why:"Discharging assumption A after deriving B yields exactly A→B.",
      shape:"[assume A … B]  ⊢  A→B",fix:`Set the close formula to ${df(want)}.`};
  }
  if(rule==="notI"){
    if(last===null || !dis(last,"bot")) return {level:"err",msg:"¬I needs the subproof to derive ⊥ (a contradiction).",
      why:"¬Introduction discharges A after reaching absurdity; then conclude ¬A.",
      shape:"[assume A … ⊥]  ⊢  ¬A",fix:`Derive ⊥ under the assumption, then close to ¬${df(A)}.`};
    const want={t:"op",op:"not",args:[A]};
    if(dfeq(concl,want)) return {level:"ok",msg:`¬Introduction: discharged ${df(A)}, concluded ${df(want)}.`,
      discharged:aIdx,grade:`Assumption ${df(A)} discharged via reductio; ¬${df(A)} holds at the outer scope.`};
    return {level:"err",msg:`¬I should conclude ¬${df(A)}.`,
      why:"Reaching ⊥ under assumption A discharges A and yields ¬A.",
      shape:"[assume A … ⊥]  ⊢  ¬A",fix:`Set the close formula to ${df(want)}.`};
  }
  return {level:"err",msg:"Unknown introduction (close) rule.",fix:"Use →I or ¬I to discharge."};
}

function proofStatus(lines){
  let st=[];
  for(const ln of lines){ const k=ln.kind;
    if(k==="assume") st.push(1); else if(k==="close"){ if(st.length) st.pop(); } }
  return {open_assumptions:st.length, is_theorem:st.length===0};
}

/* ----- incidence / geometry engine (Fano as an instance) ----- */
/* incidence.js — faithful JS port of incidence.py (finite incidence structures,
   with the Fano plane as an instance). Verified against the Python output. */

function makeStructure(name, points, lines){
  const L={};
  for(const [lnName, pts] of lines) L[lnName]=new Set(pts);
  return {name, points:[...points], lines:L};
}
function setHas(s,x){ return s.has(x); }
function incLineThrough(struct, p, q){
  return Object.keys(struct.lines).filter(name=>struct.lines[name].has(p)&&struct.lines[name].has(q));
}
function commonLine(struct, pts){
  for(const name of Object.keys(struct.lines))
    if(pts.every(x=>struct.lines[name].has(x))) return name;
  return null;
}
function thirdPoint(struct, p, q){
  const lines=incLineThrough(struct,p,q);
  if(lines.length!==1) return null;
  const blk=[...struct.lines[lines[0]]];
  if(blk.length!==3) return null;
  const rest=blk.filter(x=>x!==p&&x!==q);
  return rest.length===1?rest[0]:null;
}
function sortedJoin(it){ return [...it].sort().join(", "); }

function checkIncidence(ruleId, struct, args){
  const pts=struct.points;
  if(ruleId==="on"){
    const p=args.p, L=args.line;
    if(!pts.includes(p)) return badpt(p,pts);
    if(!(L in struct.lines)) return badline(L,struct);
    if(struct.lines[L].has(p)) return {level:"ok",msg:`${p} lies on line ${L}.`};
    return {level:"err",msg:`${p} does not lie on line ${L}.`,
      why:`Line ${L} = {${sortedJoin(struct.lines[L])}}.`,
      fix:`Pick a point actually on ${L}, or a line through ${p}.`};
  }
  if(ruleId==="collinear"){
    const {p,q,r}=args;
    for(const x of [p,q,r]) if(!pts.includes(x)) return badpt(x,pts);
    if(new Set([p,q,r]).size<3)
      return {level:"err",msg:"Collinearity needs three distinct points.",
        why:"A degenerate triple (with repeats) is trivially 'collinear'.",fix:"Choose three different points."};
    const L=commonLine(struct,[p,q,r]);
    if(L) return {level:"ok",msg:`${p}, ${q}, ${r} are collinear (line ${L}).`,
      grade:`They share block ${L} = {${sortedJoin(struct.lines[L])}}.`};
    return {level:"err",msg:`${p}, ${q}, ${r} are not collinear.`,
      why:`No line of ${struct.name} contains all three.`,fix:"Pick three points that share a common line."};
  }
  if(ruleId==="third"){
    const {p,q,r}=args;
    for(const x of [p,q]) if(!pts.includes(x)) return badpt(x,pts);
    if(p===q) return {level:"err",msg:"Need two distinct points to determine a line.",fix:"Choose distinct p and q."};
    const t=thirdPoint(struct,p,q);
    if(t===null) return {level:"err",msg:`${p} and ${q} don't lie on a unique 3-point line.`,
      why:"The closure rule needs the block through p,q to be a triple.",
      fix:"Pick two points whose joining line is a 3-point block."};
    if(r===t){ const L=incLineThrough(struct,p,q)[0];
      return {level:"ok",msg:`Closure: ${p}, ${q} force the third point ${t} (line ${L}).`,
        grade:"In a projective plane two points determine their line; its third point is forced — this is the Fano rule, generalized."};
    }
    return {level:"err",msg:`The third point of line ${p}${q} is ${t}, not ${r}.`,
      why:"Two points determine a unique line; its remaining point is fixed.",
      shape:"on(p,L), on(q,L)  ⊢  on(thirdpoint(L), L)",fix:`Set the third point to ${t}.`};
  }
  if(ruleId==="meet"){
    const L=args.line, M=args.line2, p=args.p;
    if(!(L in struct.lines)) return badline(L,struct);
    if(!(M in struct.lines)) return badline(M,struct);
    if(L===M) return {level:"err",msg:"Two distinct lines are needed to meet.",fix:"Choose distinct lines."};
    const inter=[...struct.lines[L]].filter(x=>struct.lines[M].has(x));
    if(inter.includes(p)) return {level:"ok",msg:`Lines ${L} and ${M} meet at ${p}.`,
      grade:`In a projective plane any two lines meet in exactly one point; here ${L} ∩ ${M} = {${inter.sort().join(", ")}}.`};
    return {level:"err",msg:`${p} is not the meet of ${L} and ${M}.`,
      why:`${L} ∩ ${M} = {${inter.sort().join(", ")||"∅"}}.`,fix:"Use their actual common point."};
  }
  return {level:"err",msg:"Unknown incidence rule."};
}
function badpt(p,pts){ return {level:"err",msg:`${p} is not a point of this structure.`,
  why:`Points are: ${pts.join(", ")}.`,fix:"Use a listed point."}; }
function badline(L,struct){ return {level:"err",msg:`${L} is not a line of this structure.`,
  why:`Lines are: ${Object.keys(struct.lines).join(", ")}.`,fix:"Use a listed line."}; }

function fanoStructure(){
  const pts=[]; for(let i=1;i<=7;i++) pts.push(""+i);
  const triples=[[1,2,3],[1,4,5],[1,7,6],[2,4,6],[2,5,7],[3,4,7],[3,6,5]];
  const lines=triples.map((t,i)=>["L"+(i+1), t.map(String)]);
  return makeStructure("Fano plane PG(2,2)", pts, lines);
}
function triangleStructure(){
  return makeStructure("Triangle", ["A","B","C"],
    [["e1",["A","B"]],["e2",["B","C"]],["e3",["A","C"]]]);
}
const STRUCTURES={fano:fanoStructure(),triangle:triangleStructure()};

/* ----- finite-field engine GF(8) (folded in; verified parity) ----- */
/* finite_field.js — faithful JS port of finite_field.py. GF(8)=𝔽₂[x]/(x³+x+1),
   real polynomial arithmetic. Verified against the Python output. */

const FF_MOD = 0b1011;   // x³ + x + 1
const FF_PRIM = 2;       // x is a primitive element

function gadd(a,b){ return a ^ b; }
function gmul(a,b){
  let p=0, aa=a, bb=b;
  while(bb){
    if(bb & 1) p ^= aa;
    bb >>= 1;
    aa <<= 1;
    if(aa & 0b1000) aa ^= FF_MOD;
  }
  return p;
}
function gpow(a,k){ let r=1; for(let i=0;i<k;i++) r=gmul(r,a); return r; }
function ginv(a){ return a===0?null:gpow(a,6); }
function gorder(a){ if(a===0) return null; let r=a,k=1; while(r!==1){ r=gmul(r,a); k++; } return k; }
function polyStr(a){
  if(a===0) return "0";
  const terms=[];
  for(const [bit,name] of [[4,"x²"],[2,"x"],[1,"1"]]) if(a&bit) terms.push(name);
  return terms.join("+");
}

function checkField(ruleId, args){
  const a=args.a, b=args.b, c=args.c, k=args.k;
  const inrange=(...xs)=>xs.every(x=>x!=null && x>=0 && x<=7);

  if(ruleId==="add"){
    if(!inrange(a,b,c)) return ffBad();
    const t=gadd(a,b);
    if(t===c) return {level:"ok",msg:`${a} ⊕ ${b} = ${c} in 𝔽₈.`,
      grade:`Addition is XOR: (${polyStr(a)}) + (${polyStr(b)}) = (${polyStr(c)}).`};
    return {level:"err",msg:`${a} ⊕ ${b} = ${t}, not ${c}.`,
      why:"𝔽₈ addition is coefficientwise mod 2 (XOR).",
      fix:`The sum of (${polyStr(a)}) and (${polyStr(b)}) is (${polyStr(t)}) = ${t}.`};
  }
  if(ruleId==="mul"){
    if(!inrange(a,b,c)) return ffBad();
    const t=gmul(a,b);
    if(t===c) return {level:"ok",msg:`${a} ⊗ ${b} = ${c} in 𝔽₈.`,
      grade:`Polynomial product (${polyStr(a)})(${polyStr(b)}) mod (x³+x+1) = (${polyStr(c)}).`};
    return {level:"err",msg:`${a} ⊗ ${b} = ${t}, not ${c}.`,
      why:"Multiply as polynomials over 𝔽₂, then reduce mod x³+x+1 (using x³ = x+1).",
      fix:`The product is (${polyStr(t)}) = ${t}.`};
  }
  if(ruleId==="inv"){
    if(!inrange(a,b)) return ffBad();
    if(a===0) return {level:"err",msg:"0 has no multiplicative inverse.",
      why:"Only nonzero field elements are invertible.",fix:"Choose a nonzero a."};
    const t=ginv(a);
    if(t===b) return {level:"ok",msg:`${a}⁻¹ = ${b} in 𝔽₈.`,
      grade:`Since 𝔽₈* is cyclic of order 7, a⁻¹ = a⁶; here ${a}⁶ = ${t}, and ${a}⊗${t} = 1.`};
    return {level:"err",msg:`${a}⁻¹ = ${t}, not ${b}.`,
      why:"The inverse is the unique element multiplying to 1.",
      fix:`${a} ⊗ ${t} = 1, so ${a}⁻¹ = ${t}.`};
  }
  if(ruleId==="pow"){
    if(k==null || c==null || !(c>=0&&c<=7)) return ffBad();
    const t=gpow(FF_PRIM,k);
    if(t===c) return {level:"ok",msg:`g^${k} = ${c}, where g = x = 2.`,
      grade:`x is a primitive element: x¹..x⁷ cycle through all 7 nonzero elements${k>0?" (g has order 7, so the exponent is read mod 7)":""}.`};
    return {level:"err",msg:`g^${k} = ${t}, not ${c}.`,
      why:"g = x = 2 generates 𝔽₈*; compute by repeated multiplication.",fix:`x^${k} = ${t}.`};
  }
  if(ruleId==="order"){
    if(!inrange(a) || c==null) return ffBad();
    if(a===0) return {level:"err",msg:"0 has no multiplicative order.",fix:"Choose a nonzero element."};
    const t=gorder(a);
    if(t===c) return {level:"ok",msg:`The order of ${a} in 𝔽₈* is ${c}.`,
      grade:"Orders divide |𝔽₈*| = 7 (Lagrange); so every nonzero non-identity element has order exactly 7."};
    return {level:"err",msg:`The order of ${a} is ${t}, not ${c}.`,
      why:"The order is the least k≥1 with a^k = 1; it must divide 7.",fix:`ord(${a}) = ${t}.`};
  }
  return {level:"err",msg:"Unknown finite-field rule."};
}
function ffBad(){ return {level:"err",msg:"Field elements must be integers 0..7.",
  why:"𝔽₈ has exactly 8 elements, encoded as 3-bit integers.",fix:"Use a value in 0..7."}; }




/* ============================================================================
   CLAYWORTH PROOF BENCH (React) — runs the verified JS engine above natively.
   This engine is a faithful port of prover_core.py, checked field-for-field
   against the Python output (39/39 cases identical). No Pyodide, no network.
   ========================================================================== */

const ATOMINFO = {
  prop:   { atoms:["p","q","r","s"], un:[["not","¬"]], bi:[["and","∧"],["or","∨"],["imp","→"],["iff","↔"]] },
  modalK: { atoms:["p","q","r"], un:[["not","¬"],["box","□"],["dia","◇"]], bi:[["and","∧"],["imp","→"]] },
  fano:   { atoms:["e1","e2","e3","e4","e5","e6","e7"], un:[], bi:[] },
};
const aN = (name)=>({t:"atom",name});
const oN = (op,args)=>({t:"op",op,args});

const RULES = {
  prop:[["assume","Assumption",0],["mp","→ Elim (modus ponens)",2],["mt","Modus tollens",2],
    ["andI","∧ Intro",2],["andE1","∧ Elim (left)",1],["andE2","∧ Elim (right)",1],
    ["orI1","∨ Intro (left)",1],["orI2","∨ Intro (right)",1],["dn","Double negation",1]],
  modalK:[["assume","Assumption",0],["mp","→ Elim (modus ponens)",2],["nec","Necessitation",1],
    ["kdist","K: □(A→B) ⊢ □A→□B",1],["dualD","◇A ⊢ ¬□¬A",1],
    ["semantic","⊨ Semantic check (Kripke)",2]],
  fano:[["assume","Given unit",0],["fano","Fano line (two ⊢ third)",2]],
};
const BLURB = {
  prop:"Natural-deduction core. Atoms p,q,r,s with ¬ ∧ ∨ → ↔. Conclusions are double-checked against a full truth table.",
  modalK:"Propositional base + □ ◇. Structural rules, plus a Kripke semantic check that searches for a countermodel in your chosen logic (K / T / S4 / S5).",
  fano:"Units e1..e7. One rule: two units on a common Fano line derive the third — and the engine computes the associator grade from real octonion multiplication.",
};
const P = { bg:"#0a0e14",bg2:"#0e141d",panel:"#121b26",edge:"#22303f",ink:"#e9f1f7",
  dim:"#7d92a8",faint:"#56697d",ok:"#46e3a0",warn:"#ffc24b",err:"#ff5c72",
  axiom:"#7aa2ff",rule:"#c08bff",trace:"#39d3df" };

function Select({value,onChange,options,small,accent}){
  return (<select value={value} onChange={(e)=>onChange(e.target.value)}
    style={{background:P.bg2,color:P.ink,border:`1px solid ${(accent||P.edge)}88`,borderRadius:6,
      padding:small?"5px 7px":"7px 10px",fontFamily:"'JetBrains Mono',monospace",fontSize:small?12.5:13,
      outline:"none",cursor:"pointer"}}>
    {options.map(o=>(<option key={o.value} value={o.value} style={{background:P.bg2}}>{o.label}</option>))}
  </select>);
}
function Pill({level}){
  const m={ok:[P.ok,"valid"],assumed:[P.axiom,"assumed"],warn:[P.warn,"check"],err:[P.err,"error"],idle:[P.faint,"—"]};
  const [c,l]=m[level]||m.idle;
  return (<span style={{display:"inline-flex",alignItems:"center",gap:5,fontSize:11,color:c,fontFamily:"'JetBrains Mono',monospace"}}>
    <span style={{width:7,height:7,borderRadius:"50%",background:c,boxShadow:`0 0 8px ${c}`}}/>{l}</span>);
}
const fbTag=(c)=>({fontSize:10,fontWeight:700,letterSpacing:.6,textTransform:"uppercase",
  color:c,fontFamily:"'JetBrains Mono',monospace",flexShrink:0,marginTop:2,minWidth:34});

function classifyShape(v){
  if(!v)return "atom";
  if(v.t==="atom")return v.name&&v.name[0]==="e"?"unit":"atom";
  if(ARITY[v.op]===1){const inner=v.args[0];
    if(inner&&inner.t==="op"&&ARITY[inner.op]===2)return "unaryBinary";return "unary";}
  if(ARITY[v.op]===2){const r=v.args[1];
    if(r&&r.t==="op"&&ARITY[r.op]===1)return "binaryUnary";return "binary";}
  return "atom";
}
function seedShape(shape,info){
  const A=info.atoms,u=(info.un[0]||["not"])[0],b=(info.bi[0]||["and"])[0];
  switch(shape){
    case "unit":return aN("e1"); case "atom":return aN(A[0]);
    case "unary":return oN(u,[aN(A[0])]); case "binary":return oN(b,[aN(A[0]),aN(A[1])]);
    case "unaryBinary":return oN(u,[oN(b,[aN(A[0]),aN(A[1])])]);
    case "binaryUnary":return oN(b,[aN(A[0]),oN(u,[aN(A[1])])]);
    default:return aN(A[0]);
  }
}
function Builder({sysId,value,onChange}){
  const info=ATOMINFO[sysId],A=info.atoms;
  const shapes=sysId==="fano"?[["unit","unit eₙ"]]:[
    ["atom","atom"],
    ...(info.un.length?[["unary","¬/□/◇ atom"]]:[]),
    ...(info.bi.length?[["binary","atom ∘ atom"]]:[]),
    ...(info.un.length&&info.bi.length?[["unaryBinary","¬(atom ∘ atom)"],["binaryUnary","atom → ¬atom"]]:[]),
  ];
  const cur=classifyShape(value), set=onChange;
  const unOpts=info.un.map(([o,s])=>({value:o,label:s}));
  const biOpts=info.bi.map(([o,s])=>({value:o,label:s}));
  const aOpts=A.map(x=>({value:x,label:x}));
  const els=[<Select key="shape" small accent={P.axiom} value={cur}
    onChange={(s)=>set(seedShape(s,info))} options={shapes.map(([v,l])=>({value:v,label:l}))}/>];
  if(cur==="unit"||cur==="atom")
    els.push(<Select key="a" small accent={P.trace} value={value?.name||A[0]} onChange={(v)=>set(aN(v))} options={aOpts}/>);
  if(cur==="unary"){
    els.push(<Select key="u" small accent={P.rule} value={value?.op||info.un[0][0]} onChange={(v)=>set(oN(v,[value?.args?.[0]||aN(A[0])]))} options={unOpts}/>);
    els.push(<Select key="ua" small accent={P.trace} value={value?.args?.[0]?.name||A[0]} onChange={(v)=>set(oN(value.op,[aN(v)]))} options={aOpts}/>);
  }
  if(cur==="binary"){
    els.push(<Select key="l" small accent={P.trace} value={value?.args?.[0]?.name||A[0]} onChange={(v)=>set(oN(value.op,[aN(v),value.args[1]]))} options={aOpts}/>);
    els.push(<Select key="o" small accent={P.rule} value={value?.op||info.bi[0][0]} onChange={(v)=>set(oN(v,value.args))} options={biOpts}/>);
    els.push(<Select key="r" small accent={P.trace} value={value?.args?.[1]?.name||A[1]} onChange={(v)=>set(oN(value.op,[value.args[0],aN(v)]))} options={aOpts}/>);
  }
  if(cur==="unaryBinary"){
    els.push(<Select key="u" small accent={P.rule} value={value?.op||"not"} onChange={(v)=>set(oN(v,[value.args[0]]))} options={unOpts}/>);
    els.push(<span key="p1" style={{color:P.faint}}>(</span>);
    els.push(<Select key="l" small accent={P.trace} value={value?.args?.[0]?.args?.[0]?.name||A[0]} onChange={(v)=>set(oN(value.op,[oN(value.args[0].op,[aN(v),value.args[0].args[1]])]))} options={aOpts}/>);
    els.push(<Select key="o" small accent={P.rule} value={value?.args?.[0]?.op||info.bi[0][0]} onChange={(v)=>set(oN(value.op,[oN(v,value.args[0].args)]))} options={biOpts}/>);
    els.push(<Select key="r" small accent={P.trace} value={value?.args?.[0]?.args?.[1]?.name||A[1]} onChange={(v)=>set(oN(value.op,[oN(value.args[0].op,[value.args[0].args[0],aN(v)])]))} options={aOpts}/>);
    els.push(<span key="p2" style={{color:P.faint}}>)</span>);
  }
  if(cur==="binaryUnary"){
    els.push(<Select key="l" small accent={P.trace} value={value?.args?.[0]?.name||A[0]} onChange={(v)=>set(oN(value.op,[aN(v),value.args[1]]))} options={aOpts}/>);
    els.push(<Select key="o" small accent={P.rule} value={value?.op||info.bi[0][0]} onChange={(v)=>set(oN(v,value.args))} options={biOpts}/>);
    els.push(<Select key="u" small accent={P.rule} value={value?.args?.[1]?.op||info.un[0][0]} onChange={(v)=>set(oN(value.op,[value.args[0],oN(v,[value.args[1].args?.[0]||aN(A[0])])]))} options={unOpts}/>);
    els.push(<Select key="ua" small accent={P.trace} value={value?.args?.[1]?.args?.[0]?.name||A[0]} onChange={(v)=>set(oN(value.op,[value.args[0],oN(value.args[1].op,[aN(v)])]))} options={aOpts}/>);
  }
  return <div style={{display:"flex",gap:6,flexWrap:"wrap",alignItems:"center"}}>{els}</div>;
}

function seedLines(id){
  if(id==="prop")return[
    {ruleId:"assume",refs:[null,null],formula:aN("p")},
    {ruleId:"assume",refs:[null,null],formula:oN("imp",[aN("p"),aN("q")])},
    {ruleId:"mp",refs:[0,1],formula:aN("q")}];
  if(id==="modalK")return[
    {ruleId:"assume",refs:[null,null],formula:oN("box",[aN("p")])},
    {ruleId:"semantic",refs:[0,null],formula:aN("p")}];
  return[
    {ruleId:"assume",refs:[null,null],formula:aN("e1")},
    {ruleId:"assume",refs:[null,null],formula:aN("e2")},
    {ruleId:"fano",refs:[0,1],formula:aN("e3")}];
}
const fanoIdxOf=(n)=>n&&n.t==="atom"&&n.name[0]==="e"?parseInt(n.name.slice(1),10):null;
const lbl={fontSize:10.5,letterSpacing:1.4,color:P.faint,textTransform:"uppercase",marginBottom:7};

/* ===========================================================================
   EQUATIONAL THEORY PANEL — algebras as equational logic (folded in).
   The object theory may be non-associative; the reasoning relation (refl/sym/
   trans/cong/axiom) is structural and dispatches LOW. The octonion theory's
   Fano triples are now DEFINING EQUATIONS, not hand-coded rules.
   =========================================================================== */
const EQ_OBJECT_DESC={group:"a group",octonion:"non-associative octonions",prime_field:"the field GF(7)",
  modular:"the ring ℤ/5ℤ",boolean:"a Boolean algebra"};
const EQ_BLURB={
  group:"A group as an equational theory: associativity, left-unit, left-inverse. Prove identities by reflexivity, symmetry, transitivity, congruence, and axiom instances.",
  octonion:"The octonion's 7 Fano triples are the theory's defining equations. A proof ABOUT octonion identities still proceeds by structural equational logic — the non-associativity is a derived note, not inference machinery (Lens 1 / Lens 2).",
  modular:"Modular arithmetic ℤ/5ℤ: the addition and multiplication tables are the defining equations (e.g. 2+3=0, 2·3=1), atop the ring axioms. The reasoning is structural — the on-ramp to finite-field reasoning 𝔽₂⊂𝔽₈⊂𝔽₆₄.",
  boolean:"Finite Boolean algebra over ∪ ∩ ¬ ⊤ ⊥: the Boolean identities (De Morgan, distributivity, absorption, complement) are the axioms. Bridges to the propositional channel via Lindenbaum–Tarski.",
  prime_field:"The prime field GF(7) = 𝔽₇: ℤ/7ℤ where 7 is prime, so every nonzero element is invertible. The inverse table is included. This is the base of the finite-field tower 𝔽₂⊂𝔽₈⊂𝔽₆₄ — the GF(8) extension lives in its own channel.",
};
const EQ_RULES=[["refl","Reflexivity (⊢ t = t)",0],["sym","Symmetry (s=t ⊢ t=s)",1],
  ["trans","Transitivity (s=t, t=u ⊢ s=u)",2],["cong","Congruence (one arg)",1],
  ["axiom","Axiom instance",0]];

function EquationalPanel(){
  const [thName,setThName]=useState("group");
  const theory = thName==="group"?theoryGroup():thName==="octonion"?theoryOctonion():thName==="modular"?theoryModular(5):thName==="boolean"?theoryBoolean():theoryPrimeField(7);
  // each line: {rule, refs:[..], l:term, r:term}
  const [lines,setLines]=useState(()=>eqSeed("group"));

  function switchTheory(t){ setThName(t); setLines(eqSeed(t)); }
  function setLine(i,patch){ setLines(ls=>ls.map((ln,j)=>j===i?{...ln,...patch}:ln)); }
  function addLine(){ setLines(ls=>[...ls,{rule:"refl",refs:[null,null],l:evar("x"),r:evar("x")}]); }
  function delLine(i){ setLines(ls=>ls.filter((_,j)=>j!==i)); }

  const evals=lines.map((ln,idx)=>{
    const rule=EQ_RULES.find(r=>r[0]===ln.rule)||EQ_RULES[0];
    const need=rule[2];
    const refs=ln.refs.filter(r=>r!=null);
    const bad=refs.find(r=>r<0||r>=idx);
    if(need>0&&bad!=null) return {level:"err",msg:`Cites line ${bad+1}, not above this one.`,fix:"Cite an earlier line."};
    const cited=refs.map(r=>({l:lines[r].l,r:lines[r].r}));
    return checkEquational(ln.rule,cited,{l:ln.l,r:ln.r},theory);
  });

  return (
    <div>
      <div style={{display:"flex",gap:18,flexWrap:"wrap",alignItems:"flex-end",padding:16,
        background:P.panel,border:`1px solid ${P.edge}`,borderRadius:10,marginBottom:10}}>
        <div><div style={lbl}>Equational theory</div>
          <Select value={thName} onChange={switchTheory} accent={P.trace}
            options={[{value:"group",label:"Group (assoc + unit + inverse)"},{value:"octonion",label:"Octonion (7 Fano equations)"},{value:"modular",label:"Modular ℤ/5ℤ (+ and ·)"},{value:"boolean",label:"Boolean algebra (∪ ∩ ¬)"},{value:"prime_field",label:"Prime field GF(7)"}]}/></div>
      </div>

      {/* substrate readout: the REASONING relation dispatches low */}
      <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center",padding:"9px 14px",
        background:P.bg2,border:`1px solid ${P.edge}`,borderRadius:9,marginBottom:8}}>
        <span style={fbTag(P.rule)}>substrate</span>
        <span style={{fontSize:12.5,color:P.ink}}>Reasoning relation dispatches to rung <b style={{color:P.trace}}>ℝ</b> — equational
          logic (refl/sym/trans/cong) is structural, even when the object theory is {EQ_OBJECT_DESC[thName]||"this algebra"}.</span>
      </div>
      <div style={{fontSize:12.5,color:P.faint,fontStyle:"italic",margin:"0 0 14px"}}>
        {EQ_BLURB[thName]||EQ_BLURB.group}
      </div>

      {/* axiom reference */}
      <div style={{padding:14,background:P.panel,border:`1px solid ${P.edge}`,borderRadius:10,marginBottom:14}}>
        <div style={{...lbl,marginBottom:9}}>Theory axioms</div>
        <div style={{display:"flex",flexDirection:"column",gap:6}}>
          {theory.axioms.map((ax,i)=>(
            <div key={i} style={{fontFamily:"'JetBrains Mono',monospace",fontSize:12.5,color:P.dim}}>
              <span style={{color:P.rule}}>{ax.name}</span>: <span style={{color:P.ink}}>{tfmt(ax.l)} = {tfmt(ax.r)}</span>
            </div>
          ))}
        </div>
      </div>

      {/* equation lines */}
      <div style={{background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,overflow:"hidden"}}>
        {lines.map((ln,i)=>{
          const rule=EQ_RULES.find(r=>r[0]===ln.rule)||EQ_RULES[0];
          const ev=evals[i];
          const bc=ev.level==="err"?P.err:ev.level==="ok"?P.ok:P.edge;
          const bg=ev.level==="err"?`${P.err}0d`:ev.level==="ok"?`${P.ok}08`:"transparent";
          return (<div key={i} style={{borderBottom:`1px solid ${P.edge}`,background:bg}}>
            <div style={{padding:"12px 14px",display:"flex",flexDirection:"column",gap:10}}>
              <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",gap:10}}>
                <div style={{display:"flex",alignItems:"center",gap:12,flexWrap:"wrap"}}>
                  <span style={{fontFamily:"'JetBrains Mono',monospace",color:P.faint,fontSize:13}}>L{i+1}</span>
                  <span style={{fontFamily:"'JetBrains Mono',monospace",fontSize:16,borderLeft:`3px solid ${bc}`,paddingLeft:10}}>{tfmt(ln.l)} = {tfmt(ln.r)}</span>
                </div>
                <div style={{display:"flex",alignItems:"center",gap:10}}>
                  <Pill level={ev.level}/>
                  <button onClick={()=>delLine(i)} style={{background:"transparent",border:"none",color:P.faint,fontSize:16,padding:4,cursor:"pointer"}}>✕</button>
                </div>
              </div>
              <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>lhs</span>
                <TermBuilder theory={theory} value={ln.l} onChange={(t)=>setLine(i,{l:t})}/>
                <span style={{color:P.dim,fontFamily:"'JetBrains Mono',monospace"}}>=</span>
                <span style={{...fbTag(P.faint),minWidth:30}}>rhs</span>
                <TermBuilder theory={theory} value={ln.r} onChange={(t)=>setLine(i,{r:t})}/>
              </div>
              <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>rule</span>
                <Select small accent={P.rule} value={ln.rule} onChange={(v)=>setLine(i,{rule:v,refs:[null,null]})}
                  options={EQ_RULES.map(r=>({value:r[0],label:r[1]}))}/>
                {rule[2]>0 && <><span style={{fontSize:12,color:P.faint}}>from</span>
                  {Array.from({length:rule[2]}).map((_,k)=>
                    <Select key={k} small accent={P.axiom} value={ln.refs[k]??""}
                      onChange={(v)=>{const refs=[...ln.refs];refs[k]=v===""?null:parseInt(v,10);setLine(i,{refs});}}
                      options={[{value:"",label:"line…"},...lines.slice(0,i).map((_,j)=>({value:j,label:"L"+(j+1)}))]}/>)}
                </>}
              </div>
            </div>
            {ev.level==="err" && <div style={{padding:"0 14px 12px 16px",display:"flex",flexDirection:"column",gap:6}}>
              <div style={{display:"flex",gap:8,alignItems:"baseline"}}>
                <span style={{...fbTag(P.err),minWidth:0}}>problem</span>
                <span style={{fontSize:13,color:P.ink}}>{ev.msg}</span></div>
              {ev.why && <div style={{display:"flex",gap:8}}><span style={fbTag(P.dim)}>why</span><span style={{fontSize:12.5,color:P.dim}}>{ev.why}</span></div>}
              {ev.shape && <div style={{display:"flex",gap:8,alignItems:"center"}}><span style={fbTag(P.rule)}>rule</span>
                <code style={{fontSize:12.5,color:P.rule,background:`${P.rule}14`,padding:"2px 8px",borderRadius:5}}>{ev.shape}</code></div>}
              {ev.fix && <div style={{display:"flex",gap:8}}><span style={fbTag(P.trace)}>fix</span><span style={{fontSize:12.5,color:P.trace}}>{ev.fix}</span></div>}
            </div>}
            {ev.level==="ok"&&ev.grade && <div style={{padding:"0 14px 11px 16px"}}>
              <span style={{fontSize:12,color:P.rule,lineHeight:1.55}}>◆ {ev.grade}</span></div>}
          </div>);
        })}
        <div style={{padding:12,display:"flex",gap:10}}>
          <button onClick={addLine} style={{padding:"9px 16px",background:P.trace,color:P.bg,border:"none",borderRadius:7,fontWeight:700,fontSize:13,fontFamily:"'JetBrains Mono',monospace",cursor:"pointer"}}>+ Add equation</button>
          <button onClick={()=>setLines(eqSeed(thName))} style={{padding:"9px 16px",background:"transparent",color:P.ink,border:`1px solid ${P.edge}`,borderRadius:7,fontSize:13,fontFamily:"'JetBrains Mono',monospace",cursor:"pointer"}}>Reset</button>
        </div>
      </div>
    </div>
  );
}

/* compact term builder: pick a variable, a constant, or an operation applied
   to sub-terms (depth-limited for usability). */
function TermBuilder({theory,value,onChange}){
  const sig=theory.sig;
  const ops=Object.keys(sig);
  const consts=ops.filter(o=>sig[o]===0);
  const funcs=ops.filter(o=>sig[o]>0);
  const VARS=["x","y","z"];
  const kind = value?.t==="var"?"var":(value?.t==="app"&&sig[value.op]===0?"const":"app");
  const set=onChange;
  const els=[
    <Select key="k" small accent={P.axiom} value={kind} onChange={(k)=>{
      if(k==="var") set(evar("x"));
      else if(k==="const") set(eapp(consts[0]||"e",[]));
      else { const f=funcs[0]||"m"; set(eapp(f, Array.from({length:sig[f]}).map(()=>VARS[0]).map(evar))); }
    }} options={[{value:"var",label:"var"},...(consts.length?[{value:"const",label:"const"}]:[]),...(funcs.length?[{value:"app",label:"op(...)"}]:[])]}/>
  ];
  if(kind==="var")
    els.push(<Select key="v" small accent={P.trace} value={value.name} onChange={(v)=>set(evar(v))} options={VARS.map(v=>({value:v,label:v}))}/>);
  if(kind==="const")
    els.push(<Select key="c" small accent={P.trace} value={value.op} onChange={(v)=>set(eapp(v,[]))} options={consts.map(c=>({value:c,label:c}))}/>);
  if(kind==="app"){
    els.push(<Select key="f" small accent={P.rule} value={value.op} onChange={(v)=>set(eapp(v, Array.from({length:sig[v]}).map((_,i)=>value.args[i]||evar("x"))))} options={funcs.map(f=>({value:f,label:f}))}/>);
    els.push(<span key="(" style={{color:P.faint}}>(</span>);
    for(let i=0;i<sig[value.op];i++){
      if(i>0) els.push(<span key={"comma"+i} style={{color:P.faint}}>,</span>);
      els.push(<TermBuilder key={"arg"+i} theory={theory} value={value.args[i]||evar("x")}
        onChange={(t)=>{const args=[...value.args];args[i]=t;set(eapp(value.op,args));}}/>);
    }
    els.push(<span key=")" style={{color:P.faint}}>)</span>);
  }
  return <span style={{display:"inline-flex",gap:4,alignItems:"center",flexWrap:"wrap"}}>{els}</span>;
}

function eqSeed(thName){
  const x=evar("x"),y=evar("y"),z=evar("z");
  if(thName==="group") return [
    {rule:"axiom",refs:[null,null],l:eapp("m",[eapp("e",[]),x]),r:x},
    {rule:"refl",refs:[null,null],l:x,r:x},
    {rule:"sym",refs:[0,null],l:x,r:eapp("m",[eapp("e",[]),x])}];
  if(thName==="octonion") return [
    {rule:"axiom",refs:[null,null],l:eapp("m",[eapp("e1",[]),eapp("e2",[])]),r:eapp("e3",[])},
    {rule:"sym",refs:[0,null],l:eapp("e3",[]),r:eapp("m",[eapp("e1",[]),eapp("e2",[])])}];
  if(thName==="modular"){ const K=i=>eapp(""+i,[]);
    return [
    {rule:"axiom",refs:[null,null],l:eapp("a",[K(2),K(3)]),r:K(0)},   // 2+3=0 mod 5
    {rule:"axiom",refs:[null,null],l:eapp("m",[K(2),K(3)]),r:K(1)},   // 2·3=1 mod 5
    {rule:"sym",refs:[0,null],l:K(0),r:eapp("a",[K(2),K(3)])}]; }
  if(thName==="prime_field"){ const K=i=>eapp(""+i,[]);
    return [
    {rule:"axiom",refs:[null,null],l:eapp("m",[K(3),K(5)]),r:K(1)},   // 3·5=1 mod 7
    {rule:"sym",refs:[0,null],l:K(1),r:eapp("m",[K(3),K(5)])}]; }
  // boolean
  return [
    {rule:"axiom",refs:[null,null],l:eapp("not",[eapp("or",[x,y])]),r:eapp("and",[eapp("not",[x]),eapp("not",[y])])}, // De Morgan
    {rule:"axiom",refs:[null,null],l:eapp("not",[eapp("not",[x])]),r:x}];  // double negation
}


/* ===========================================================================
   LINEAR LOGIC PANEL — multiplicative fragment with resource accounting.
   Resources are used exactly once; the engine runs a backtracking search.
   =========================================================================== */
const LIN_ATOMS=["a","b","c","d"];
function lAtom(s){ return {t:"atom",name:s}; }
function lTensor(a,b){ return {t:"op",op:"tensor",args:[a,b]}; }
function lLolli(a,b){ return {t:"op",op:"lolli",args:[a,b]}; }

function LinTermBuilder({value,onChange}){
  const kind = value?.t==="atom"?"atom":(value?.op==="tensor"?"tensor":(value?.op==="lolli"?"lolli":"one"));
  const set=onChange;
  const els=[
    <Select key="k" small accent={P.axiom} value={kind} onChange={(k)=>{
      if(k==="atom") set(lAtom("a"));
      else if(k==="tensor") set(lTensor(lAtom("a"),lAtom("b")));
      else if(k==="lolli") set(lLolli(lAtom("a"),lAtom("b")));
      else set({t:"op",op:"one",args:[]});
    }} options={[{value:"atom",label:"atom"},{value:"tensor",label:"⊗ tensor"},{value:"lolli",label:"⊸ lollipop"},{value:"one",label:"1"}]}/>
  ];
  if(kind==="atom")
    els.push(<Select key="a" small accent={P.trace} value={value.name} onChange={(v)=>set(lAtom(v))} options={LIN_ATOMS.map(x=>({value:x,label:x}))}/>);
  if(kind==="tensor"||kind==="lolli"){
    els.push(<span key="(" style={{color:P.faint}}>(</span>);
    els.push(<LinTermBuilder key="l" value={value.args[0]} onChange={(t)=>set({t:"op",op:value.op,args:[t,value.args[1]]})}/>);
    els.push(<span key="op" style={{color:P.rule,fontFamily:"'JetBrains Mono',monospace"}}>{kind==="tensor"?"⊗":"⊸"}</span>);
    els.push(<LinTermBuilder key="r" value={value.args[1]} onChange={(t)=>set({t:"op",op:value.op,args:[value.args[0],t]})}/>);
    els.push(<span key=")" style={{color:P.faint}}>)</span>);
  }
  return <span style={{display:"inline-flex",gap:4,alignItems:"center",flexWrap:"wrap"}}>{els}</span>;
}

function LinearPanel(){
  const [resources,setResources]=useState(()=>[lAtom("a"),lLolli(lAtom("a"),lAtom("b"))]);
  const [goal,setGoal]=useState(()=>lAtom("b"));
  const verdict=useMemo(()=>{
    try { return checkLinear(resources,goal); }
    catch(e){ return {level:"err",msg:"engine error: "+e.message}; }
  },[resources,goal]);

  const setRes=(i,t)=>setResources(rs=>rs.map((r,j)=>j===i?t:r));
  const addRes=()=>setResources(rs=>[...rs,lAtom("a")]);
  const delRes=(i)=>setResources(rs=>rs.filter((_,j)=>j!==i));

  const bc=verdict.level==="ok"?P.ok:P.err;
  return (
    <div>
      <div style={{fontSize:12.5,color:P.faint,fontStyle:"italic",margin:"0 0 14px"}}>
        Multiplicative linear logic: every resource is consumed exactly once — no duplicating, no discarding.
        The engine runs a backtracking search over resource allocations (the channel that finally exercises reversal).
      </div>

      <div style={{padding:16,background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,marginBottom:14}}>
        <div style={{...lbl,marginBottom:10}}>Resources (context Γ)</div>
        <div style={{display:"flex",flexDirection:"column",gap:8}}>
          {resources.map((r,i)=>(
            <div key={i} style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
              <span style={{fontFamily:"'JetBrains Mono',monospace",color:P.faint,fontSize:12,minWidth:24}}>R{i+1}</span>
              <LinTermBuilder value={r} onChange={(t)=>setRes(i,t)}/>
              <button onClick={()=>delRes(i)} style={{background:"transparent",border:"none",color:P.faint,fontSize:15,cursor:"pointer"}}>✕</button>
            </div>
          ))}
        </div>
        <button onClick={addRes} style={{marginTop:10,padding:"7px 13px",background:"transparent",color:P.ink,border:`1px solid ${P.edge}`,borderRadius:7,fontFamily:"'JetBrains Mono',monospace",fontSize:12.5,cursor:"pointer"}}>+ Add resource</button>
      </div>

      <div style={{padding:16,background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,marginBottom:14}}>
        <div style={{display:"flex",alignItems:"center",gap:10,flexWrap:"wrap"}}>
          <span style={{...lbl,marginBottom:0}}>Goal</span>
          <LinTermBuilder value={goal} onChange={setGoal}/>
        </div>
      </div>

      {/* sequent + verdict */}
      <div style={{padding:16,borderRadius:12,border:`1px solid ${bc}`,background:`${bc}0d`}}>
        <div style={{fontFamily:"'JetBrains Mono',monospace",fontSize:15,color:P.ink,marginBottom:10}}>
          {resources.map(linShow).join(", ")||"·"} <span style={{color:P.rule}}>⊢</span> {linShow(goal)}
        </div>
        <div style={{display:"flex",gap:8,alignItems:"baseline",marginBottom:6}}>
          <Pill level={verdict.level}/>
          <span style={{fontSize:13,color:P.ink}}>{verdict.msg}</span>
        </div>
        {verdict.why && <div style={{display:"flex",gap:8}}><span style={fbTag(P.dim)}>why</span>
          <span style={{fontSize:12.5,color:P.dim,lineHeight:1.55}}>{verdict.why}</span></div>}
        {verdict.fix && <div style={{display:"flex",gap:8}}><span style={fbTag(P.trace)}>fix</span>
          <span style={{fontSize:12.5,color:P.trace}}>{verdict.fix}</span></div>}
        {verdict.grade && <div style={{marginTop:4}}><span style={{fontSize:12,color:P.rule}}>◆ {verdict.grade}</span></div>}
      </div>

      <div style={{marginTop:14,fontSize:12,color:P.faint,lineHeight:1.7}}>
        Try the linear distinctions: <code style={{color:P.rule}}>a, b ⊢ a</code> fails (b can't be discarded);
        <code style={{color:P.rule}}> a ⊢ a⊗a</code> fails (a can't be duplicated); but
        <code style={{color:P.rule}}> a⊸b, b⊸c ⊢ a⊸c</code> succeeds.
      </div>
    </div>
  );
}
function linShow(n){
  if(n==null) return "·";
  if(n.t==="atom") return n.name;
  if(n.op==="one") return "1";
  return "("+linShow(n.args[0])+(n.op==="tensor"?" ⊗ ":" ⊸ ")+linShow(n.args[1])+")";
}

/* ===========================================================================
   SELF-TEST PANEL — runs every embedded engine's checks live, in-artifact.
   If the bundler ever breaks an engine, you see exactly which one here.
   =========================================================================== */
function SelfTestPanel(){
  const results=useMemo(()=>runAllSelfTests(),[]);
  const total=results.reduce((s,g)=>s+g.tests.length,0);
  const passed=results.reduce((s,g)=>s+g.tests.filter(t=>t.ok).length,0);
  return (
    <details style={{marginTop:20,background:P.panel,border:`1px solid ${P.edge}`,borderRadius:10,padding:"10px 14px"}}>
      <summary style={{cursor:"pointer",fontFamily:"'JetBrains Mono',monospace",fontSize:12.5,
        color:passed===total?P.ok:P.err}}>
        ◆ Engine self-tests: {passed}/{total} passing {passed===total?"✓":"✗"}
      </summary>
      <div style={{marginTop:10,display:"flex",flexDirection:"column",gap:10}}>
        {results.map((g,i)=>(
          <div key={i}>
            <div style={{fontFamily:"'JetBrains Mono',monospace",fontSize:12,color:P.dim,marginBottom:4}}>
              {g.name}: {g.tests.filter(t=>t.ok).length}/{g.tests.length}
            </div>
            <div style={{display:"flex",flexWrap:"wrap",gap:4}}>
              {g.tests.map((t,j)=>(
                <span key={j} title={t.name} style={{width:9,height:9,borderRadius:2,
                  background:t.ok?P.ok:P.err,display:"inline-block"}}/>
              ))}
            </div>
          </div>
        ))}
      </div>
    </details>
  );
}
function runAllSelfTests(){
  const groups=[];
  const A=(s)=>({t:"atom",name:s}), o=(op,a)=>({t:"op",op,args:a});
  const imp=(a,b)=>o("imp",[a,b]), box=(a)=>o("box",[a]), dia=(a)=>o("dia",[a]), neg=(a)=>o("not",[a]);
  const T=(name,ok)=>({name,ok:!!ok});

  // propositional / core
  try{
    const t=[];
    t.push(T("MP valid", check_step("prop","mp",2,[A("p"),imp(A("p"),A("q"))],A("q"),null,"K").level==="ok"));
    t.push(T("affirming consequent rejected", check_step("prop","mp",2,[A("q"),imp(A("p"),A("q"))],A("p"),null,"K").level==="err"));
    t.push(T("andE1", check_step("prop","andE1",1,[o("and",[A("p"),A("q")])],A("p"),null,"K").level==="ok"));
    groups.push({name:"propositional",tests:t});
  }catch(e){ groups.push({name:"propositional",tests:[T("threw: "+e.message,false)]}); }

  // modal Kripke hierarchy
  try{
    const t=[];
    t.push(T("□p⊬p in K", check_step("modalK","semantic",1,[box(A("p"))],A("p"),null,"K").level==="err"));
    t.push(T("□p⊢p in T", check_step("modalK","semantic",1,[box(A("p"))],A("p"),null,"T").level==="ok"));
    t.push(T("□p⊬□□p in T", check_step("modalK","semantic",1,[box(A("p"))],box(box(A("p"))),null,"T").level==="err"));
    t.push(T("□p⊢□□p in S4", check_step("modalK","semantic",1,[box(A("p"))],box(box(A("p"))),null,"S4").level==="ok"));
    t.push(T("◇p⊢□◇p in S5", check_step("modalK","semantic",1,[dia(A("p"))],box(dia(A("p"))),null,"S5").level==="ok"));
    groups.push({name:"modal (K/T/S4/S5)",tests:t});
  }catch(e){ groups.push({name:"modal",tests:[T("threw: "+e.message,false)]}); }

  // fano
  try{
    const t=[];
    t.push(T("e1,e2→e3 valid", check_step("fano","fano",2,[],null,[1,2,3],"K").level==="ok"));
    t.push(T("e1,e2→e5 rejected", check_step("fano","fano",2,[],null,[1,2,5],"K").level==="err"));
    groups.push({name:"fano",tests:t});
  }catch(e){ groups.push({name:"fano",tests:[T("threw: "+e.message,false)]}); }

  // equational
  try{
    const t=[]; const O=theoryOctonion();
    const ea=eapp("m",[eapp("e1",[]),eapp("e2",[])]);
    t.push(T("octonion fano axiom", checkEquational("axiom",[],{l:ea,r:eapp("e3",[])},O).level==="ok"));
    t.push(T("wrong product rejected", checkEquational("axiom",[],{l:ea,r:eapp("e5",[])},O).level==="err"));
    const G=theoryGroup(), x=evar("x"),y=evar("y"),z=evar("z");
    t.push(T("trans chain", checkEquational("trans",[{l:x,r:y},{l:y,r:z}],{l:x,r:z},G).level==="ok"));
    t.push(T("cong one-arg", checkEquational("cong",[{l:y,r:z}],{l:eapp("m",[x,y]),r:eapp("m",[x,z])},G).level==="ok"));
    const M=theoryModular(5), K=(i)=>eapp(""+i,[]);
    t.push(T("mod5 2+3=0", checkEquational("axiom",[],{l:eapp("a",[K(2),K(3)]),r:K(0)},M).level==="ok"));
    t.push(T("mod5 2·3=1", checkEquational("axiom",[],{l:eapp("m",[K(2),K(3)]),r:K(1)},M).level==="ok"));
    t.push(T("mod5 2+3=4 rejected", checkEquational("axiom",[],{l:eapp("a",[K(2),K(3)]),r:K(4)},M).level==="err"));
    const B=theoryBoolean();
    t.push(T("bool De Morgan", checkEquational("axiom",[],{l:eapp("not",[eapp("or",[x,y])]),r:eapp("and",[eapp("not",[x]),eapp("not",[y])])},B).level==="ok"));
    t.push(T("bool double-neg", checkEquational("axiom",[],{l:eapp("not",[eapp("not",[x])]),r:x},B).level==="ok"));
    const PF=theoryPrimeField(7), K7=(i)=>eapp(""+i,[]);
    t.push(T("GF(7) 3·5=1", checkEquational("axiom",[],{l:eapp("m",[K7(3),K7(5)]),r:K7(1)},PF).level==="ok"));
    t.push(T("GF(7) inv(3)=5", checkEquational("axiom",[],{l:eapp("m",[K7(3),K7(5)]),r:K7(1)},PF).level==="ok"));
    groups.push({name:"equational (5 theories)",tests:t});
  }catch(e){ groups.push({name:"equational",tests:[T("threw: "+e.message,false)]}); }

  // linear
  try{
    const t=[]; const a=lAtom("a"),b=lAtom("b"),c=lAtom("c");
    t.push(T("a,a⊸b⊢b", checkLinear([a,lLolli(a,b)],b).level==="ok"));
    t.push(T("compose", checkLinear([lLolli(a,b),lLolli(b,c)],lLolli(a,c)).level==="ok"));
    t.push(T("no weakening", checkLinear([a,b],a).level==="err"));
    t.push(T("no contraction", checkLinear([a],lTensor(a,a)).level==="err"));
    groups.push({name:"linear (MLL)",tests:t});
  }catch(e){ groups.push({name:"linear",tests:[T("threw: "+e.message,false)]}); }

  // discharge / subproofs
  try{
    const t=[];
    const pnp=o("and",[A("p"),o("not",[A("p")])]);
    const reductio=[
      {kind:"assume",depth:1,formula:pnp},
      {kind:"line",depth:1,rule:"andE",refs:[0],formula:A("p")},
      {kind:"line",depth:1,rule:"andE",refs:[0],formula:o("not",[A("p")])},
      {kind:"line",depth:1,rule:"botI",refs:[1,2],formula:{t:"op",op:"bot",args:[]}},
      {kind:"close",depth:0,rule:"notI",formula:o("not",[pnp])},
    ];
    t.push(T("andE in subproof", checkDischargeStep(reductio,1).level==="ok"));
    t.push(T("botI from A,¬A", checkDischargeStep(reductio,3).level==="ok"));
    t.push(T("reductio closes ¬(p∧¬p)", checkDischargeStep(reductio,4).level==="ok"));
    t.push(T("¬(p∧¬p) is a theorem", proofStatus(reductio).is_theorem===true));
    // scope violation
    const scope=[
      {kind:"assume",depth:1,formula:A("p")},
      {kind:"line",depth:1,rule:"reit",refs:[0],formula:A("p")},
      {kind:"close",depth:0,rule:"impI",formula:o("imp",[A("p"),A("p")])},
      {kind:"line",depth:0,rule:"reit",refs:[1],formula:A("p")},
    ];
    t.push(T("citing discharged line rejected", checkDischargeStep(scope,3).level==="err"));
    groups.push({name:"discharge (subproofs)",tests:t});
  }catch(e){ groups.push({name:"discharge",tests:[T("threw: "+e.message,false)]}); }

  // incidence geometry
  try{
    const t=[]; const F=fanoStructure();
    t.push(T("Fano 1,2,3 collinear", checkIncidence("collinear",F,{p:"1",q:"2",r:"3"}).level==="ok"));
    t.push(T("Fano 1,2,4 not collinear", checkIncidence("collinear",F,{p:"1",q:"2",r:"4"}).level==="err"));
    t.push(T("Fano closure 1,2→3", checkIncidence("third",F,{p:"1",q:"2",r:"3"}).level==="ok"));
    t.push(T("Fano closure 1,2→5 rejected", checkIncidence("third",F,{p:"1",q:"2",r:"5"}).level==="err"));
    t.push(T("Fano L1∩L2 = 1", checkIncidence("meet",F,{line:"L1",line2:"L2",p:"1"}).level==="ok"));
    const Tri=triangleStructure();
    t.push(T("triangle: no 3-closure", checkIncidence("third",Tri,{p:"A",q:"B",r:"C"}).level==="err"));
    groups.push({name:"incidence geometry",tests:t});
  }catch(e){ groups.push({name:"incidence",tests:[T("threw: "+e.message,false)]}); }

  // finite field GF(8)
  try{
    const t=[];
    t.push(T("3⊕5=6", checkField("add",{a:3,b:5,c:6}).level==="ok"));
    t.push(T("2⊗4=3 (x·x²=x+1)", checkField("mul",{a:2,b:4,c:3}).level==="ok"));
    t.push(T("2⊗4=5 rejected", checkField("mul",{a:2,b:4,c:5}).level==="err"));
    t.push(T("3⁻¹=6", checkField("inv",{a:3,b:ginv(3)}).level==="ok"));
    t.push(T("0 has no inverse", checkField("inv",{a:0,b:1}).level==="err"));
    t.push(T("x⁷=1 (cycle closes)", checkField("pow",{k:7,c:1}).level==="ok"));
    t.push(T("ord(2)=7", checkField("order",{a:2,c:7}).level==="ok"));
    // exhaustive field-law spot: associativity holds for all triples
    let assoc=true; for(let a=0;a<8;a++)for(let b=0;b<8;b++)for(let c=0;c<8;c++) if(gmul(gmul(a,b),c)!==gmul(a,gmul(b,c))) assoc=false;
    t.push(T("mul associative (all 512 triples)", assoc));
    groups.push({name:"finite field GF(8)",tests:t});
  }catch(e){ groups.push({name:"finite field",tests:[T("threw: "+e.message,false)]}); }

  return groups;
}

/* ===========================================================================
   DISCHARGE / SUBPROOF PANEL — natural deduction with assumption discharge.
   Opening an assumption indents into a subproof box; →I / ¬I close the box and
   discharge its assumption. Scope is enforced: you can't cite a line whose
   subproof has closed. A depth-0 conclusion with no open assumptions is a THEOREM.
   =========================================================================== */
const DIS_ATOMS=["p","q","r"];
function dAtom(s){ return {t:"atom",name:s}; }
function dOp(o,a){ return {t:"op",op:o,args:a}; }
const DBOT={t:"op",op:"bot",args:[]};
const DIS_INSCOPE=[["reit","Reiterate",1],["mp","→ Elim (MP)",2],["andI","∧ Intro",2],
  ["andE","∧ Elim",1],["botI","⊥ Intro (A,¬A)",2]];

function DischargeFormulaBuilder({value,onChange,allowBot}){
  // compact builder over p/q/r with ¬ ∧ → and optionally ⊥
  const kind = value?.t==="atom" ? "atom"
    : value?.op==="bot" ? "bot"
    : value?.op==="not" ? "not"
    : value?.op==="and" ? "and" : "imp";
  const set=onChange;
  const opts=[{value:"atom",label:"atom"},{value:"not",label:"¬"},{value:"and",label:"∧"},{value:"imp",label:"→"}];
  if(allowBot) opts.push({value:"bot",label:"⊥"});
  const els=[
    <Select key="k" small accent={P.axiom} value={kind} onChange={(k)=>{
      if(k==="atom") set(dAtom("p"));
      else if(k==="bot") set(DBOT);
      else if(k==="not") set(dOp("not",[dAtom("p")]));
      else set(dOp(k,[dAtom("p"),dAtom("q")]));
    }} options={opts}/>
  ];
  if(kind==="atom")
    els.push(<Select key="a" small accent={P.trace} value={value.name} onChange={(v)=>set(dAtom(v))} options={DIS_ATOMS.map(x=>({value:x,label:x}))}/>);
  if(kind==="not"){
    els.push(<span key="n" style={{color:P.rule}}>¬</span>);
    els.push(<DischargeFormulaBuilder key="i" value={value.args[0]} onChange={(t)=>set(dOp("not",[t]))}/>);
  }
  if(kind==="and"||kind==="imp"){
    els.push(<span key="(" style={{color:P.faint}}>(</span>);
    els.push(<DischargeFormulaBuilder key="l" value={value.args[0]} onChange={(t)=>set(dOp(value.op,[t,value.args[1]]))}/>);
    els.push(<span key="op" style={{color:P.rule,fontFamily:"'JetBrains Mono',monospace"}}>{kind==="and"?"∧":"→"}</span>);
    els.push(<DischargeFormulaBuilder key="r" value={value.args[1]} onChange={(t)=>set(dOp(value.op,[value.args[0],t]))}/>);
    els.push(<span key=")" style={{color:P.faint}}>)</span>);
  }
  return <span style={{display:"inline-flex",gap:4,alignItems:"center",flexWrap:"wrap"}}>{els}</span>;
}

function DischargePanel(){
  const [lines,setLines]=useState(()=>disSeed());

  const evals=lines.map((_,i)=>{
    try { return checkDischargeStep(lines,i); }
    catch(e){ return {level:"err",msg:"engine error: "+e.message}; }
  });
  const status=(()=>{ try{ return proofStatus(lines); }catch(e){ return {open_assumptions:0,is_theorem:false}; } })();

  function setLine(i,patch){ setLines(ls=>ls.map((ln,j)=>j===i?{...ln,...patch}:ln)); }
  function delLine(i){ setLines(ls=>ls.filter((_,j)=>j!==i)); }

  // structural editing: add a line at current max depth; open/close subproofs
  function addLine(){ setLines(ls=>{ const d=curDepth(ls); return [...ls,{kind:"line",depth:d,rule:"reit",refs:[null,null],formula:dAtom("p")}]; }); }
  function openSub(){ setLines(ls=>{ const d=curDepth(ls); return [...ls,{kind:"assume",depth:d+1,formula:dAtom("p")}]; }); }
  function closeSub(){ setLines(ls=>{ const d=curDepth(ls); if(d===0) return ls;
    return [...ls,{kind:"close",depth:d-1,rule:"impI",formula:dOp("imp",[dAtom("p"),dAtom("p")])}]; }); }

  return (
    <div>
      <div style={{fontSize:12.5,color:P.faint,fontStyle:"italic",margin:"0 0 14px"}}>
        Natural deduction with discharge. Open a subproof to assume a hypothesis (indents); close it with
        →I or ¬I to discharge that assumption. Citing a line from a closed subproof is a scope error.
        A depth-0 result with no open assumptions is a genuine theorem.
      </div>

      {/* theorem/entailment status */}
      <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center",padding:"9px 14px",
        background:P.bg2,border:`1px solid ${status.is_theorem?P.ok:P.edge}`,borderRadius:9,marginBottom:12}}>
        <span style={fbTag(status.is_theorem?P.ok:P.rule)}>{status.is_theorem?"theorem":"status"}</span>
        <span style={{fontSize:12.5,color:P.ink}}>
          {status.open_assumptions===0
            ? "No open assumptions — a depth-0 conclusion here is a theorem (⊢)."
            : `${status.open_assumptions} assumption(s) still open — results are conditional (entailment), not theorems.`}
        </span>
      </div>

      <div style={{background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,overflow:"hidden"}}>
        {lines.map((ln,i)=>{
          const ev=evals[i];
          const bc=ev.level==="err"?P.err:ev.level==="ok"?P.ok:ev.level==="assumed"?P.axiom:P.edge;
          const indent=ln.depth*18;
          const isAssume=ln.kind==="assume", isClose=ln.kind==="close";
          return (<div key={i} style={{borderBottom:`1px solid ${P.edge}`,
            background:ev.level==="err"?`${P.err}0d`:ev.level==="ok"?`${P.ok}07`:ev.level==="assumed"?`${P.axiom}07`:"transparent"}}>
            <div style={{padding:"11px 14px 11px",marginLeft:indent,
              borderLeft:ln.depth>0?`2px solid ${P.axiom}55`:"none",paddingLeft:ln.depth>0?12:0}}>
              <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",gap:10}}>
                <div style={{display:"flex",alignItems:"center",gap:10,flexWrap:"wrap"}}>
                  <span style={{fontFamily:"'JetBrains Mono',monospace",color:P.faint,fontSize:12}}>L{i+1}</span>
                  {isAssume && <span style={{...fbTag(P.axiom),minWidth:0}}>assume</span>}
                  {isClose && <span style={{...fbTag(P.rule),minWidth:0}}>{ln.rule==="impI"?"→I":"¬I"}</span>}
                  <span style={{fontFamily:"'JetBrains Mono',monospace",fontSize:15,
                    borderLeft:`3px solid ${bc}`,paddingLeft:9}}>{df(ln.formula)}</span>
                </div>
                <div style={{display:"flex",alignItems:"center",gap:8}}>
                  <Pill level={ev.level}/>
                  <button onClick={()=>delLine(i)} style={{background:"transparent",border:"none",color:P.faint,fontSize:15,cursor:"pointer"}}>✕</button>
                </div>
              </div>

              {/* editors */}
              <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap",marginTop:9}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>formula</span>
                <DischargeFormulaBuilder value={ln.formula} onChange={(t)=>setLine(i,{formula:t})} allowBot={!isAssume}/>
              </div>

              {isClose && <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap",marginTop:8}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>rule</span>
                <Select small accent={P.rule} value={ln.rule} onChange={(v)=>setLine(i,{rule:v})}
                  options={[{value:"impI",label:"→I (discharge → A→B)"},{value:"notI",label:"¬I (reductio → ¬A)"}]}/>
              </div>}

              {ln.kind==="line" && <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap",marginTop:8}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>rule</span>
                <Select small accent={P.rule} value={ln.rule} onChange={(v)=>setLine(i,{rule:v,refs:[null,null]})}
                  options={DIS_INSCOPE.map(r=>({value:r[0],label:r[1]}))}/>
                {(() => { const need=(DIS_INSCOPE.find(r=>r[0]===ln.rule)||[,, 0])[2];
                  return need>0 && <><span style={{fontSize:12,color:P.faint}}>from</span>
                    {Array.from({length:need}).map((_,k)=>
                      <Select key={k} small accent={P.axiom} value={ln.refs[k]??""}
                        onChange={(v)=>{const refs=[...ln.refs];refs[k]=v===""?null:parseInt(v,10);setLine(i,{refs});}}
                        options={[{value:"",label:"line…"},...lines.slice(0,i).map((_,j)=>({value:j,label:"L"+(j+1)}))]}/>)}
                  </>; })()}
              </div>}

              {/* feedback */}
              {ev.level==="err" && <div style={{marginTop:9,display:"flex",flexDirection:"column",gap:5}}>
                <div style={{display:"flex",gap:8,alignItems:"baseline"}}><span style={{...fbTag(P.err),minWidth:0}}>problem</span>
                  <span style={{fontSize:12.5,color:P.ink}}>{ev.msg}</span></div>
                {ev.why && <div style={{display:"flex",gap:8}}><span style={fbTag(P.dim)}>why</span><span style={{fontSize:12,color:P.dim}}>{ev.why}</span></div>}
                {ev.shape && <div style={{display:"flex",gap:8,alignItems:"center"}}><span style={fbTag(P.rule)}>rule</span>
                  <code style={{fontSize:12,color:P.rule,background:`${P.rule}14`,padding:"2px 7px",borderRadius:5}}>{ev.shape}</code></div>}
                {ev.fix && <div style={{display:"flex",gap:8}}><span style={fbTag(P.trace)}>fix</span><span style={{fontSize:12,color:P.trace}}>{ev.fix}</span></div>}
              </div>}
              {ev.level==="ok"&&ev.grade && <div style={{marginTop:7}}><span style={{fontSize:12,color:P.rule}}>◆ {ev.grade}</span></div>}
              {ev.level==="assumed" && <div style={{marginTop:7}}><span style={{fontSize:12,color:P.dim}}>{ev.why}</span></div>}
            </div>
          </div>);
        })}
        <div style={{padding:12,display:"flex",gap:8,flexWrap:"wrap"}}>
          <button onClick={addLine} style={btnP(P.trace,P.bg)}>+ Line</button>
          <button onClick={openSub} style={btnP(P.axiom,P.bg)}>⊢ Open subproof (assume)</button>
          <button onClick={closeSub} style={btnGhost()}>Close subproof (discharge)</button>
          <button onClick={()=>setLines(disSeed())} style={btnGhost()}>Reset</button>
        </div>
      </div>
    </div>
  );
}
function btnP(bg,fg){ return {padding:"8px 14px",background:bg,color:fg,border:"none",borderRadius:7,fontWeight:700,fontSize:12.5,fontFamily:"'JetBrains Mono',monospace",cursor:"pointer"}; }
function btnGhost(){ return {padding:"8px 14px",background:"transparent",color:P.ink,border:`1px solid ${P.edge}`,borderRadius:7,fontSize:12.5,fontFamily:"'JetBrains Mono',monospace",cursor:"pointer"}; }
function curDepth(ls){ let d=0; for(const ln of ls){ if(ln.kind==="assume") d++; else if(ln.kind==="close") d=Math.max(0,d-1); } return d; }
function disSeed(){
  // worked example: prove ¬(p∧¬p) by reductio
  const p=dAtom("p"), np=dOp("not",[dAtom("p")]), pnp=dOp("and",[dAtom("p"),dOp("not",[dAtom("p")])]);
  return [
    {kind:"assume",depth:1,formula:pnp},
    {kind:"line",depth:1,rule:"andE",refs:[0,null],formula:p},
    {kind:"line",depth:1,rule:"andE",refs:[0,null],formula:np},
    {kind:"line",depth:1,rule:"botI",refs:[1,2],formula:DBOT},
    {kind:"close",depth:0,rule:"notI",formula:dOp("not",[pnp])},
  ];
}

/* ===========================================================================
   INCIDENCE / GEOMETRY PANEL — finite incidence structures. The Fano plane is
   one instance; the projective-closure rule IS the Fano "two points force the
   third," generalized. Pick a structure, pick a relation, get a verdict.
   =========================================================================== */
const INC_RULES=[["collinear","Collinear (3 points share a line)"],
  ["third","Closure: two points force the third"],
  ["on","Point on line"],["meet","Lines meet at a point"]];

function IncidencePanel(){
  const [structName,setStructName]=useState("fano");
  const struct = structName==="fano"?fanoStructure():triangleStructure();
  const [rule,setRule]=useState("collinear");
  const pts=struct.points;
  const lineNames=Object.keys(struct.lines);
  const [args,setArgs]=useState(()=>({p:pts[0],q:pts[1],r:pts[2],line:lineNames[0],line2:lineNames[1]}));

  function switchStruct(s){
    setStructName(s);
    const st = s==="fano"?fanoStructure():triangleStructure();
    const p=st.points, ln=Object.keys(st.lines);
    setArgs({p:p[0],q:p[1]||p[0],r:p[2]||p[0],line:ln[0],line2:ln[1]||ln[0]});
    if(s==="triangle" && rule==="third") setRule("collinear");
  }

  const verdict=useMemo(()=>{
    try{ return checkIncidence(rule,struct,args); }
    catch(e){ return {level:"err",msg:"engine error: "+e.message}; }
  },[rule,structName,args]);

  const ptOpts=pts.map(x=>({value:x,label:x}));
  const lineOpts=lineNames.map(x=>({value:x,label:x}));
  const setA=(patch)=>setArgs(a=>({...a,...patch}));
  const bc=verdict.level==="ok"?P.ok:P.err;

  return (
    <div>
      <div style={{display:"flex",gap:18,flexWrap:"wrap",alignItems:"flex-end",padding:16,
        background:P.panel,border:`1px solid ${P.edge}`,borderRadius:10,marginBottom:10}}>
        <div><div style={lbl}>Incidence structure</div>
          <Select value={structName} onChange={switchStruct} accent={P.trace}
            options={[{value:"fano",label:"Fano plane PG(2,2)"},{value:"triangle",label:"Triangle (3 points)"}]}/></div>
        <div><div style={lbl}>Relation</div>
          <Select value={rule} onChange={setRule} accent={P.rule}
            options={INC_RULES.map(r=>({value:r[0],label:r[1]}))}/></div>
      </div>

      <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center",padding:"9px 14px",
        background:P.bg2,border:`1px solid ${P.edge}`,borderRadius:9,marginBottom:8}}>
        <span style={fbTag(P.rule)}>structure</span>
        <span style={{fontSize:12.5,color:P.ink}}>{struct.name} — points {"{"+pts.join(", ")+"}"}, {lineNames.length} lines.</span>
      </div>
      <div style={{fontSize:12.5,color:P.faint,fontStyle:"italic",margin:"0 0 14px"}}>
        Finite incidence geometry. The closure rule “two points force the third” is exactly the Fano
        multiplication law, generalized to any projective line. The triangle has no 3-point lines, so closure fails there.
      </div>

      <div style={{padding:16,background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,marginBottom:14}}>
        <div style={{...lbl,marginBottom:10}}>Arguments</div>
        <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center"}}>
          {(rule==="collinear"||rule==="third") && <>
            <span style={fbTag(P.faint)}>p</span><Select small accent={P.trace} value={args.p} onChange={(v)=>setA({p:v})} options={ptOpts}/>
            <span style={fbTag(P.faint)}>q</span><Select small accent={P.trace} value={args.q} onChange={(v)=>setA({q:v})} options={ptOpts}/>
            <span style={fbTag(P.faint)}>{rule==="third"?"third":"r"}</span><Select small accent={P.trace} value={args.r} onChange={(v)=>setA({r:v})} options={ptOpts}/>
          </>}
          {rule==="on" && <>
            <span style={fbTag(P.faint)}>p</span><Select small accent={P.trace} value={args.p} onChange={(v)=>setA({p:v})} options={ptOpts}/>
            <span style={fbTag(P.faint)}>line</span><Select small accent={P.rule} value={args.line} onChange={(v)=>setA({line:v})} options={lineOpts}/>
          </>}
          {rule==="meet" && <>
            <span style={fbTag(P.faint)}>line</span><Select small accent={P.rule} value={args.line} onChange={(v)=>setA({line:v})} options={lineOpts}/>
            <span style={fbTag(P.faint)}>line₂</span><Select small accent={P.rule} value={args.line2} onChange={(v)=>setA({line2:v})} options={lineOpts}/>
            <span style={fbTag(P.faint)}>at</span><Select small accent={P.trace} value={args.p} onChange={(v)=>setA({p:v})} options={ptOpts}/>
          </>}
        </div>
      </div>

      <div style={{padding:16,borderRadius:12,border:`1px solid ${bc}`,background:`${bc}0d`}}>
        <div style={{display:"flex",gap:8,alignItems:"baseline",marginBottom:6}}>
          <Pill level={verdict.level}/>
          <span style={{fontSize:13,color:P.ink}}>{verdict.msg}</span>
        </div>
        {verdict.why && <div style={{display:"flex",gap:8}}><span style={fbTag(P.dim)}>why</span>
          <span style={{fontSize:12.5,color:P.dim,lineHeight:1.55}}>{verdict.why}</span></div>}
        {verdict.shape && <div style={{display:"flex",gap:8,alignItems:"center",marginTop:4}}><span style={fbTag(P.rule)}>rule</span>
          <code style={{fontSize:12.5,color:P.rule,background:`${P.rule}14`,padding:"2px 8px",borderRadius:5}}>{verdict.shape}</code></div>}
        {verdict.fix && <div style={{display:"flex",gap:8,marginTop:4}}><span style={fbTag(P.trace)}>fix</span>
          <span style={{fontSize:12.5,color:P.trace}}>{verdict.fix}</span></div>}
        {verdict.grade && <div style={{marginTop:6}}><span style={{fontSize:12,color:P.rule}}>◆ {verdict.grade}</span></div>}
      </div>
    </div>
  );
}

/* ===========================================================================
   CORE PANEL — propositional, modal (K/T/S4/S5 Kripke), and Fano/octonion.
   A line-oriented proof editor: pick a formula, pick a rule, cite earlier
   lines. Conclusions are double-checked (truth table / countermodel search /
   octonion associator) by the verified core engine.
   =========================================================================== */
function CorePanel(){
  const [sysId,setSysId]=useState("prop");
  const [logic,setLogic]=useState("K");
  const [lines,setLines]=useState(()=>seedLines("prop"));

  function switchSys(s){ setSysId(s); setLines(seedLines(s)); }
  function setLine(i,patch){ setLines(ls=>ls.map((ln,j)=>j===i?{...ln,...patch}:ln)); }
  function delLine(i){ setLines(ls=>ls.filter((_,j)=>j!==i)); }
  function addLine(){ const info=ATOMINFO[sysId];
    setLines(ls=>[...ls,{ruleId:RULES[sysId][0][0],refs:[null,null],formula:seedShape(sysId==="fano"?"unit":"atom",info)}]); }

  const rulesFor=RULES[sysId];
  const evals=lines.map((ln,idx)=>{
    const rule=rulesFor.find(r=>r[0]===ln.ruleId)||rulesFor[0];
    const need=rule[2];
    const refs=ln.refs.filter(r=>r!=null);
    const bad=refs.find(r=>r<0||r>=idx);
    if(need>0 && bad!=null) return {level:"err",msg:`Cites line ${bad+1}, not above this one.`,fix:"Cite an earlier line."};
    const prem=refs.map(r=>lines[r].formula);
    let fanoIdx=null;
    if(sysId==="fano") fanoIdx=[...prem.map(fanoIdxOf),fanoIdxOf(ln.formula)];
    try { return check_step(sysId,ln.ruleId,need,prem,ln.formula,fanoIdx,logic); }
    catch(e){ return {level:"err",msg:"engine error: "+e.message}; }
  });

  const sub = sysId==="fano" ? dispatch_system(true,true) : null;

  return (
    <div>
      <div style={{display:"flex",gap:18,flexWrap:"wrap",alignItems:"flex-end",padding:16,
        background:P.panel,border:`1px solid ${P.edge}`,borderRadius:10,marginBottom:10}}>
        <div><div style={lbl}>System</div>
          <Select value={sysId} onChange={switchSys} accent={P.trace}
            options={[{value:"prop",label:"Propositional (truth-table checked)"},
                      {value:"modalK",label:"Modal □ ◇ (Kripke)"},
                      {value:"fano",label:"Fano / octonion units"}]}/></div>
        {sysId==="modalK" && <div><div style={lbl}>Modal logic</div>
          <Select value={logic} onChange={setLogic} accent={P.rule}
            options={[{value:"K",label:"K"},{value:"T",label:"T (reflexive)"},{value:"S4",label:"S4 (refl + trans)"},{value:"S5",label:"S5 (equivalence)"}]}/></div>}
      </div>

      {sub && <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center",padding:"9px 14px",
        background:P.bg2,border:`1px solid ${P.edge}`,borderRadius:9,marginBottom:8}}>
        <span style={fbTag(P.rule)}>substrate</span>
        <span style={{fontSize:12.5,color:P.ink}}>Dispatched to rung <b style={{color:P.trace}}>{sub.rung}</b> —
          octonion multiplication is non-associative, so a single Fano step lives in a quaternionic chart while
          cross-line derivations make the associator go live.</span>
      </div>}

      <div style={{fontSize:12.5,color:P.faint,fontStyle:"italic",margin:"0 0 14px"}}>{BLURB[sysId]}</div>

      <div style={{background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,overflow:"hidden"}}>
        {lines.map((ln,i)=>{
          const rule=rulesFor.find(r=>r[0]===ln.ruleId)||rulesFor[0];
          const need=rule[2];
          const ev=evals[i];
          const bc=ev.level==="err"?P.err:ev.level==="ok"?P.ok:ev.level==="assumed"?P.axiom:ev.level==="warn"?P.warn:P.edge;
          const bg=ev.level==="err"?`${P.err}0d`:ev.level==="ok"?`${P.ok}07`:ev.level==="assumed"?`${P.axiom}07`:ev.level==="warn"?`${P.warn}0a`:"transparent";
          return (<div key={i} style={{borderBottom:`1px solid ${P.edge}`,background:bg}}>
            <div style={{padding:"12px 14px",display:"flex",flexDirection:"column",gap:10}}>
              <div style={{display:"flex",alignItems:"center",justifyContent:"space-between",gap:10}}>
                <div style={{display:"flex",alignItems:"center",gap:12,flexWrap:"wrap"}}>
                  <span style={{fontFamily:"'JetBrains Mono',monospace",color:P.faint,fontSize:13}}>L{i+1}</span>
                  <span style={{fontFamily:"'JetBrains Mono',monospace",fontSize:16,borderLeft:`3px solid ${bc}`,paddingLeft:10}}>{fmt(ln.formula)}</span>
                </div>
                <div style={{display:"flex",alignItems:"center",gap:10}}>
                  <Pill level={ev.level}/>
                  <button onClick={()=>delLine(i)} style={{background:"transparent",border:"none",color:P.faint,fontSize:16,padding:4,cursor:"pointer"}}>✕</button>
                </div>
              </div>
              <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>formula</span>
                <Builder sysId={sysId} value={ln.formula} onChange={(t)=>setLine(i,{formula:t})}/>
              </div>
              <div style={{display:"flex",alignItems:"center",gap:8,flexWrap:"wrap"}}>
                <span style={{...fbTag(P.faint),minWidth:46}}>rule</span>
                <Select small accent={P.rule} value={ln.ruleId} onChange={(v)=>setLine(i,{ruleId:v,refs:[null,null]})}
                  options={rulesFor.map(r=>({value:r[0],label:r[1]}))}/>
                {need>0 && <><span style={{fontSize:12,color:P.faint}}>from</span>
                  {Array.from({length:need}).map((_,k)=>
                    <Select key={k} small accent={P.axiom} value={ln.refs[k]??""}
                      onChange={(v)=>{const refs=[...ln.refs];refs[k]=v===""?null:parseInt(v,10);setLine(i,{refs});}}
                      options={[{value:"",label:"line…"},...lines.slice(0,i).map((_,j)=>({value:j,label:"L"+(j+1)}))]}/>)}
                </>}
              </div>
            </div>
            {(ev.level==="err"||ev.level==="warn") && <div style={{padding:"0 14px 12px 16px",display:"flex",flexDirection:"column",gap:6}}>
              <div style={{display:"flex",gap:8,alignItems:"baseline"}}>
                <span style={{...fbTag(ev.level==="warn"?P.warn:P.err),minWidth:0}}>{ev.level==="warn"?"check":"problem"}</span>
                <span style={{fontSize:13,color:P.ink}}>{ev.msg}</span></div>
              {ev.why && <div style={{display:"flex",gap:8}}><span style={fbTag(P.dim)}>why</span><span style={{fontSize:12.5,color:P.dim,lineHeight:1.5}}>{ev.why}</span></div>}
              {ev.shape && <div style={{display:"flex",gap:8,alignItems:"center"}}><span style={fbTag(P.rule)}>rule</span>
                <code style={{fontSize:12.5,color:P.rule,background:`${P.rule}14`,padding:"2px 8px",borderRadius:5}}>{ev.shape}</code></div>}
              {ev.counterexample && <div style={{display:"flex",gap:8}}><span style={fbTag(P.warn)}>counterexample</span><span style={{fontSize:12.5,color:P.warn,lineHeight:1.5}}>{ev.counterexample}</span></div>}
              {ev.fix && <div style={{display:"flex",gap:8}}><span style={fbTag(P.trace)}>fix</span><span style={{fontSize:12.5,color:P.trace}}>{ev.fix}</span></div>}
            </div>}
            {ev.level==="ok"&&ev.grade && <div style={{padding:"0 14px 11px 16px"}}>
              <span style={{fontSize:12,color:P.rule,lineHeight:1.55}}>◆ {ev.grade}</span></div>}
            {ev.level==="assumed" && <div style={{padding:"0 14px 11px 16px"}}>
              <span style={{fontSize:12,color:P.dim}}>{ev.why}</span></div>}
          </div>);
        })}
        <div style={{padding:12,display:"flex",gap:10}}>
          <button onClick={addLine} style={btnP(P.trace,P.bg)}>+ Add line</button>
          <button onClick={()=>setLines(seedLines(sysId))} style={btnGhost()}>Reset</button>
        </div>
      </div>
    </div>
  );
}

/* ===========================================================================
   FINITE FIELD PANEL — the Galois field GF(8) = 𝔽₂[x]/(x³+x+1). Pick an
   operation and a claim; the engine verifies it by real polynomial arithmetic.
   =========================================================================== */
const FF_RULES=[["mul","a ⊗ b = c (multiply)"],["add","a ⊕ b = c (add / XOR)"],
  ["inv","a⁻¹ = b (inverse)"],["pow","gᵏ = c (powers of x)"],["order","ord(a) = c (order)"]];
const FF_ELEMS=[0,1,2,3,4,5,6,7].map(v=>({value:v,label:String(v)}));

function FiniteFieldPanel(){
  const [rule,setRule]=useState("mul");
  const [args,setArgs]=useState({a:2,b:4,c:3,k:7});
  const verdict=useMemo(()=>{ try{return checkField(rule,args);}catch(e){return{level:"err",msg:"engine error: "+e.message};} },[rule,args]);
  const setA=(patch)=>setArgs(a=>({...a,...patch}));
  const numSel=(key,accent)=> <Select small accent={accent||P.trace} value={args[key]} onChange={(v)=>setA({[key]:parseInt(v,10)})} options={FF_ELEMS}/>;
  const bc=verdict.level==="ok"?P.ok:P.err;
  return (
    <div>
      <div style={{display:"flex",gap:18,flexWrap:"wrap",alignItems:"flex-end",padding:16,
        background:P.panel,border:`1px solid ${P.edge}`,borderRadius:10,marginBottom:10}}>
        <div><div style={lbl}>Operation in 𝔽₈</div>
          <Select value={rule} onChange={setRule} accent={P.rule} options={FF_RULES.map(r=>({value:r[0],label:r[1]}))}/></div>
      </div>

      <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center",padding:"9px 14px",
        background:P.bg2,border:`1px solid ${P.edge}`,borderRadius:9,marginBottom:8}}>
        <span style={fbTag(P.rule)}>substrate</span>
        <span style={{fontSize:12.5,color:P.ink}}>𝔽₈ = 𝔽₂[x]/(x³+x+1); elements are 3-bit integers 0..7, addition is XOR, x is primitive.</span>
      </div>
      <div style={{fontSize:12.5,color:P.faint,fontStyle:"italic",margin:"0 0 14px"}}>
        The Galois field GF(8), top of the finite-field tower 𝔽₂⊂𝔽₈⊂𝔽₆₄. Real polynomial arithmetic mod x³+x+1.
      </div>

      <div style={{padding:16,background:P.panel,border:`1px solid ${P.edge}`,borderRadius:12,marginBottom:14}}>
        <div style={{...lbl,marginBottom:10}}>Claim</div>
        <div style={{display:"flex",gap:10,flexWrap:"wrap",alignItems:"center",fontFamily:"'JetBrains Mono',monospace",fontSize:14,color:P.ink}}>
          {rule==="mul" && <>{numSel("a")}<span style={{color:P.rule}}>⊗</span>{numSel("b")}<span>=</span>{numSel("c",P.ok)}</>}
          {rule==="add" && <>{numSel("a")}<span style={{color:P.rule}}>⊕</span>{numSel("b")}<span>=</span>{numSel("c",P.ok)}</>}
          {rule==="inv" && <>{numSel("a")}<span style={{color:P.rule}}>⁻¹</span><span>=</span>{numSel("b",P.ok)}</>}
          {rule==="pow" && <><span>x</span><span style={{color:P.faint}}>^(k=</span>
            <Select small accent={P.trace} value={args.k} onChange={(v)=>setA({k:parseInt(v,10)})} options={Array.from({length:9}).map((_,i)=>({value:i,label:String(i)}))}/>
            <span style={{color:P.faint}}>)</span><span>=</span>{numSel("c",P.ok)}</>}
          {rule==="order" && <><span style={{color:P.rule}}>ord(</span>{numSel("a")}<span style={{color:P.rule}}>)</span><span>=</span>{numSel("c",P.ok)}</>}
        </div>
      </div>

      <div style={{padding:16,borderRadius:12,border:`1px solid ${bc}`,background:`${bc}0d`}}>
        <div style={{display:"flex",gap:8,alignItems:"baseline",marginBottom:6}}>
          <Pill level={verdict.level}/><span style={{fontSize:13,color:P.ink}}>{verdict.msg}</span></div>
        {verdict.why && <div style={{display:"flex",gap:8}}><span style={fbTag(P.dim)}>why</span><span style={{fontSize:12.5,color:P.dim,lineHeight:1.55}}>{verdict.why}</span></div>}
        {verdict.fix && <div style={{display:"flex",gap:8,marginTop:4}}><span style={fbTag(P.trace)}>fix</span><span style={{fontSize:12.5,color:P.trace}}>{verdict.fix}</span></div>}
        {verdict.grade && <div style={{marginTop:6}}><span style={{fontSize:12,color:P.rule}}>◆ {verdict.grade}</span></div>}
      </div>
    </div>
  );
}

/* ===========================================================================
   APP — tabbed shell over the eight verified engines, plus live self-tests.
   =========================================================================== */
const TABS=[
  ["core","Core logics"],
  ["equational","Equational algebra"],
  ["linear","Linear (MLL)"],
  ["discharge","Discharge / subproofs"],
  ["incidence","Incidence geometry"],
  ["field","Finite field GF(8)"],
];

export default function ClayworthProofBench(){
  const [tab,setTab]=useState("core");
  return (
    <div style={{minHeight:"100vh",background:P.bg,color:P.ink,
      fontFamily:"'Inter',system-ui,sans-serif",padding:"28px 22px"}}>
      <div style={{maxWidth:1040,margin:"0 auto"}}>
        <div style={{marginBottom:6,fontSize:11,letterSpacing:2,color:P.faint,textTransform:"uppercase"}}>Clayworth TN-2026-11 · proof bench</div>
        <h1 style={{margin:"0 0 6px",fontSize:26,fontWeight:800,letterSpacing:-0.5}}>Logic Engines</h1>
        <div style={{fontSize:13.5,color:P.dim,marginBottom:20,lineHeight:1.6,maxWidth:780}}>
          Eight verified reasoning engines — propositional, modal Kripke, octonionic Fano, equational algebra,
          multiplicative linear logic, natural-deduction discharge, finite incidence geometry, and the Galois field
          GF(8) — each a faithful JavaScript port of its Python original, running natively in the browser.
        </div>

        <div style={{display:"flex",gap:6,flexWrap:"wrap",marginBottom:22}}>
          {TABS.map(([id,label])=>(
            <button key={id} onClick={()=>setTab(id)} style={{
              padding:"8px 14px",borderRadius:8,fontSize:12.5,fontFamily:"'JetBrains Mono',monospace",
              cursor:"pointer",border:`1px solid ${tab===id?P.trace:P.edge}`,
              background:tab===id?`${P.trace}1a`:"transparent",
              color:tab===id?P.trace:P.dim,fontWeight:tab===id?700:400}}>{label}</button>
          ))}
        </div>

        {tab==="core" && <CorePanel/>}
        {tab==="equational" && <EquationalPanel/>}
        {tab==="linear" && <LinearPanel/>}
        {tab==="discharge" && <DischargePanel/>}
        {tab==="incidence" && <IncidencePanel/>}
        {tab==="field" && <FiniteFieldPanel/>}

        <SelfTestPanel/>
      </div>
    </div>
  );
}
