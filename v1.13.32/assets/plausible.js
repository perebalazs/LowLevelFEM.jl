// Plausible analytics loader
(function() {
  var d = document, s = d.createElement('script');
  s.src = "https://plausible.io/js/pa-juJDmOsc0aSv-gHV3eE7C.js";
  s.async = true;
  d.head.appendChild(s);

  window.plausible = window.plausible || function() {
    (plausible.q = plausible.q || []).push(arguments);
  };
  plausible.init = plausible.init || function(i) {
    plausible.o = i || {};
  };

  plausible.init();
})();

