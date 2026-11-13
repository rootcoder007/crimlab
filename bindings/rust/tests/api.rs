// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#[test]
fn every_function_is_there() {
    let f = crimlab::jfns();
    for n in ["vodtst", "orrisk", "sntbnd", "kinidx", "selftest", "proof"] {
        assert!(f.contains(&format!("\"{n}\":")), "{n}");
    }
}

#[test]
fn results() {
    let r = crimlab::orrisk(r#"{"or": 2, "rate": 0.1, "share": 0.3}"#).unwrap();
    assert!(r.starts_with(r#"{"exposed":0.14752049189933314,"#), "{r}");
    let r = crimlab::jrun("dyadrt", r#"{"count": [[640, 160], [160, 40]], "pop": {"a": 800, "b": 200}}"#).unwrap();
    assert!(r.starts_with(r#"{"off":["a","b","a","b"],"vic":["a","a","b","b"]"#), "{r}");
    assert!(crimlab::selftest("{}").unwrap().contains(r#""nfail":0,"report":"""#));
}

#[test]
fn errors() {
    let e = crimlab::orrisk(r#"{"or": 2, "rate": 1, "share": 0.3}"#).unwrap_err();
    assert_eq!(e.0, "orrisk: `rate` must be one number strictly between 0 and 1.");
    assert!(crimlab::jrun("nope", "{}").unwrap_err().0.contains("no function nope"));
    assert!(crimlab::losbnd("[1").is_err());
    assert!(crimlab::jrun("cllr", "{\"same\":\0}").is_err());
}

#[test]
fn large_result() {
    // longer than the first buffer: the call is made again with room
    let v: Vec<String> = (1..20000).map(|i| (i as f64 / 7.0).to_string()).collect();
    let a = format!(r#"{{"same": [{}], "diff": [{}]}}"#, v.join(","), v.join(","));
    let r = crimlab::lrcal(&a).unwrap();
    assert!(r.len() > 1 << 16);
}
